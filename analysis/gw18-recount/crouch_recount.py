#!/usr/bin/env python3
"""SRA tag-preserving conversion, integer MEX packaging, and author-data subset.

No third-party Python dependencies. SRA SEQUENCE is read once per exported row;
the original primary/secondary alignment multiplicity cannot duplicate a read.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import math
import re
import sys
from pathlib import Path
from urllib.parse import quote


def opened(path, mode="rt"):
    path = Path(path)
    return gzip.open(path, mode) if path.suffix == ".gz" else path.open(mode)


def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def integer_vector(value):
    value = value.strip().strip("[]")
    if not re.fullmatch(r"\d+(?:\s*,\s*\d+)*", value):
        raise ValueError("Unexpected numeric vector: " + repr(value))
    return [int(v) for v in value.split(",")]


def sam_record(row):
    # vdb-dump -I -f csv -C READ,SPOT_GROUP,LINKAGE_GROUP,READ_LEN,READ_START,
    #                         (INSDC:quality:text:phred_33)QUALITY
    if len(row) != 7:
        raise ValueError("Expected row ID plus six fields; got %d fields" % len(row))
    rid, seq, group, linkage, lengths, starts, qual = row
    if not rid.isdigit() or int(rid) < 1:
        raise ValueError("Invalid SEQUENCE row ID: " + repr(rid))
    if integer_vector(lengths) != [len(seq)] or integer_vector(starts) != [0]:
        raise ValueError("Only the audited single-biological-read layout is supported; row " + rid)
    if not seq or not re.fullmatch("[ACGTNacgtn]+", seq):
        raise ValueError("Invalid read sequence at row " + rid)
    if len(qual) != len(seq) or any(not 33 <= ord(q) <= 126 for q in qual):
        raise ValueError("Original read qualities failed validation at row " + rid)
    fields = group.split(":")
    if len(fields) != 5 or not fields[0] or not fields[1] or not fields[2].isdigit() or not fields[4].isdigit():
        raise ValueError("Unrecognized read-group structure; explicit library mapping needed: " + repr(group))
    # The last two fields are flowcell and lane, not separate GEM libraries.
    library = ":".join(fields[:3])
    tags = {}
    for token in linkage.split("|"):
        if not token:
            continue
        key, sep, value = token.partition(":")
        if sep:
            if key in tags:
                raise ValueError("Repeated linkage tag at row " + rid)
            tags[key] = value
    if not tags.get("CB"):
        return None, "missing_CB", group, library
    if not tags.get("UB"):
        return None, "missing_UB", group, library
    cb, ub = tags["CB"], tags["UB"]
    if not re.fullmatch(r"[ACGT]{16}-[1-9]\d*", cb):
        raise ValueError("Unexpected corrected CB at row %s: %r" % (rid, cb))
    if not re.fullmatch(r"[ACGT]{10}", ub):
        if re.fullmatch(r"[ACGTN]{10}", ub):
            return None, "UMI_with_N", group, library
        raise ValueError("Expected the audited 10-base UB at row %s: %r" % (rid, ub))
    cell = quote(library, safe="") + "|" + cb
    # A real read and its original quality; unmapped SAM, with corrected labels.
    sam = "\t".join(["SRR17878454." + rid, "4", "*", "0", "0", "*", "*", "0", "0",
                     seq.upper(), qual, "CB:Z:" + cell, "UB:Z:" + ub]) + "\n"
    return sam, "tagged", group, library


def convert(args):
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    stats, groups, lengths = collections.Counter(), collections.Counter(), collections.Counter()
    minimum, maximum, previous = None, None, 0
    preview = (output / "input_preview.csv").open("w", newline="")
    writer = csv.writer(preview)
    writer.writerow(["row_id", "READ", "SPOT_GROUP", "LINKAGE_GROUP", "READ_LEN", "READ_START", "QUALITY_phred33"])
    handles, inputs = {}, {}
    try:
        for row in csv.reader(sys.stdin):
            if stats["rows"] < 50:
                writer.writerow(row)
            record, status, group, library = sam_record(row)
            rid = int(row[0])
            if rid <= previous:
                raise ValueError("SEQUENCE row IDs are repeated or out of order")
            previous = maximum = rid
            minimum = rid if minimum is None else minimum
            stats["rows"] += 1
            stats[status] += 1
            groups[(group, library, status)] += 1
            lengths[len(row[1])] += 1
            if record:
                fields = record.rstrip("\n").split("\t")
                full_cb = fields[11].split("|", 1)[1]
                sequence_cb, gem = full_cb.rsplit("-", 1)
                key = (library, gem)
                if key not in handles:
                    if len(handles) >= 64:
                        raise ValueError("Unexpectedly many GEM libraries; inspect grouping")
                    name = "lib_" + hashlib.sha256((library + "|" + gem).encode()).hexdigest()[:16]
                    handle = (output / (name + ".sam")).open("w", buffering=1024 * 1024)
                    handle.write("@HD\tVN:1.6\tSO:unsorted\n")
                    handles[key] = handle
                    inputs[key] = [name, library, gem, 0]
                # Separate SAM per physical library/GEM. STAR sees the original
                # 16-base corrected barcode; no synthetic barcode is introduced.
                fields[11] = "CB:Z:" + sequence_cb
                handles[key].write("\t".join(fields) + "\n")
                inputs[key][3] += 1
            if stats["rows"] % 1000000 == 0:
                print("Converted {:,} SEQUENCE rows; {:,} with CB/UB".format(stats["rows"], stats["tagged"]), file=sys.stderr, flush=True)
    finally:
        preview.close()
        for handle in handles.values():
            handle.close()
    with (output / "library_inputs.tsv").open("w") as handle:
        handle.write("file_id\tlibrary\tgem\treads\n")
        for value in sorted(inputs.values()):
            handle.write("%s\t%s\t%s\t%d\n" % tuple(value))
    report = dict(stats)
    report.update(first_row=minimum, last_row=maximum, read_lengths=dict(lengths),
                  tagged_fraction=stats["tagged"] / max(1, stats["rows"]))
    save_json(output / "conversion.json", report)
    with (output / "read_groups.tsv").open("w") as handle:
        handle.write("spot_group\tlibrary\tstatus\trows\n")
        for (group, library, status), n in sorted(groups.items()):
            handle.write("%s\t%s\t%s\t%d\n" % (group, library, status, n))
    if args.expected and stats["rows"] != args.expected:
        raise ValueError("Unexpected number of exported rows: %d instead of %d" % (stats["rows"], args.expected))
    if not stats["tagged"]:
        raise ValueError("No reads with valid corrected CB and UB")
    print(json.dumps(report, indent=2), flush=True)


def mex_header(handle):
    header = handle.readline().strip()
    if header.lower() not in ("%%matrixmarket matrix coordinate integer general", "%%matrixmarket matrix coordinate real general"):
        raise ValueError("Unsupported Matrix Market header: " + header)
    for line in handle:
        if line.strip() and not line.startswith("%"):
            dims = [int(v) for v in line.split()]
            if len(dims) != 3:
                raise ValueError("Invalid matrix dimensions")
            return dims
    raise ValueError("Missing matrix dimensions")


def package(args):
    source, dest = Path(args.source), Path(args.output)
    dest.mkdir(parents=True, exist_ok=True)
    features = (source / "features.tsv").read_text().splitlines()
    barcodes = (source / "barcodes.tsv").read_text().splitlines()
    if len(set(barcodes)) != len(barcodes):
        raise ValueError("Duplicate library-aware cell identifiers")
    with (source / "matrix.mtx").open() as handle:
        nr, nc, nz = mex_header(handle)
        if nr != len(features) or nc != len(barcodes):
            raise ValueError("MEX dimensions do not match features/barcodes")
        totals, genes = [0] * nc, [0] * nc
        seen = 0
        last = None
        for line in handle:
            if not line.strip():
                continue
            i, j, value = line.split()
            i, j, value = int(i), int(j), float(value)
            if not (1 <= i <= nr and 1 <= j <= nc and math.isfinite(value) and value > 0 and value.is_integer()):
                raise ValueError("Non-positive/non-integer/invalid entry in STARsolo unique-count matrix")
            # STAR's MEX is column-major. Also catches duplicate coordinates.
            key = (j, i)
            if last is not None and key <= last:
                raise ValueError("Unexpected duplicate/unsorted STARsolo matrix coordinate")
            last = key
            totals[j - 1] += int(value)
            genes[j - 1] += 1
            seen += 1
    if seen != nz:
        raise ValueError("Matrix nonzero count mismatch")
    for filename in ("matrix.mtx", "features.tsv", "barcodes.tsv"):
        with (source / filename).open("rb") as inp, gzip.open(dest / (filename + ".gz"), "wb", compresslevel=2) as out:
            while True:
                block = inp.read(1024 * 1024)
                if not block:
                    break
                out.write(block)
    from urllib.parse import unquote
    with (dest / "cell_metadata.tsv").open("w") as handle:
        handle.write("cell_id\tsample\trun\tlibrary\tCB\tnCount_RNA\tnFeature_RNA\n")
        for cell, total, detected in zip(barcodes, totals, genes):
            library, sep, cb = cell.rpartition("|")
            if not sep or not re.fullmatch(r"[ACGT]{16}-[1-9]\d*", cb):
                raise ValueError("STAR changed a library-aware cell identifier: " + cell)
            handle.write("%s\tGW18\tSRR17878454\t%s\t%s\t%d\t%d\n" % (cell, unquote(library), cb, total, detected))
    save_json(dest / "matrix_summary.json", dict(genes=nr, barcodes=nc, nonzero=nz, UMIs=sum(totals),
               nonzero_barcodes=sum(t > 0 for t in totals), filtered=False, normalized=False))


def merge(args):
    """Merge independent GEM-library matrices without combining barcode identities."""
    manifest, star_root, output = Path(args.manifest), Path(args.star_root), Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    with manifest.open() as handle:
        libraries = list(csv.DictReader(handle, delimiter="\t"))
    bodies, features, barcodes, total_nz = [], None, [], 0
    for lib in libraries:
        source = star_root / lib["file_id"] / "Solo.out" / args.feature / "raw"
        current_features = (source / "features.tsv").read_text()
        if features is None:
            features = current_features
        elif features != current_features:
            raise ValueError("Feature rows differ between GEM libraries")
        cb_rows = (source / "barcodes.tsv").read_text().splitlines()
        offset = len(barcodes)
        for cb in cb_rows:
            # STARsolo writes sequence barcodes without a Cell Ranger GEM suffix.
            if not re.fullmatch("[ACGT]{16}", cb):
                raise ValueError("Unexpected STARsolo barcode format: " + cb)
            barcodes.append(quote(lib["library"], safe="") + "|" + cb + "-" + lib["gem"])
        with (source / "matrix.mtx").open() as handle:
            nr, nc, nz = mex_header(handle)
        if nr != len(features.splitlines()) or nc != len(cb_rows):
            raise ValueError("Per-library MEX dimensions disagree")
        total_nz += nz
        bodies.append((source / "matrix.mtx", offset))
    if features is None:
        raise ValueError("No GEM-library matrices to merge")
    (output / "features.tsv").write_text(features)
    (output / "barcodes.tsv").write_text("\n".join(barcodes) + "\n")
    with (output / "matrix.mtx").open("w") as out:
        out.write("%%MatrixMarket matrix coordinate integer general\n")
        out.write("%d %d %d\n" % (len(features.splitlines()), len(barcodes), total_nz))
        for path, offset in bodies:
            with path.open() as handle:
                mex_header(handle)
                for line in handle:
                    if line.strip():
                        i, j, value = line.split()
                        out.write("%s %d %s\n" % (i, int(j) + offset, value))


def author_subset(args):
    """Stream the author's large MEX; retain GW18 columns without normalizing."""
    root, out = Path(args.directory), Path(args.output)
    out.mkdir(parents=True, exist_ok=False)
    with (root / "meta.tsv").open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"cellId", "age"}.issubset(reader.fieldnames):
            raise ValueError("Expected meta.tsv columns cellId and age")
        names = reader.fieldnames
        rows = [r for r in reader if r["age"].strip().upper() == "GW18"]
    if not rows:
        raise ValueError("No age=GW18 rows; inspect metadata labels")
    selected = {}
    for row in rows:
        cell = row["cellId"]
        if cell in selected:
            raise ValueError("Duplicate metadata cellId")
        selected[cell] = row
    with opened(root / "barcodes.tsv.gz") as handle:
        all_barcodes = [line.rstrip("\r\n").split("\t")[0] for line in handle]
    if len(set(all_barcodes)) != len(all_barcodes):
        raise ValueError("Duplicate author matrix barcodes")
    missing = set(selected).difference(all_barcodes)
    if missing:
        raise ValueError("Metadata and matrix cell IDs do not match exactly; do not guess/strip prefixes. Examples: " + repr(sorted(missing)[:3]))
    column_map, kept = {}, []
    for old, cell in enumerate(all_barcodes, 1):
        if cell in selected:
            kept.append(cell)
            column_map[old] = len(kept)
    body = out / ".entries.tmp"
    kept_nz = seen = 0
    with opened(root / "matrix.mtx.gz") as handle, body.open("w") as temp:
        nr, nc, nz = mex_header(handle)
        if nc != len(all_barcodes):
            raise ValueError("Author barcode count does not match matrix dimensions")
        for line in handle:
            if not line.strip():
                continue
            i, j, v = line.split()
            old = int(j)
            if not (1 <= int(i) <= nr and 1 <= old <= nc and math.isfinite(float(v))):
                raise ValueError("Invalid author matrix entry")
            seen += 1
            if old in column_map:
                temp.write("%s %d %s\n" % (i, column_map[old], v))
                kept_nz += 1
            if seen % 10000000 == 0:
                print("Scanned {:,}/{:,} author entries".format(seen, nz), flush=True)
    if seen != nz:
        raise ValueError("Truncated author matrix")
    with gzip.open(out / "matrix.mtx.gz", "wt", compresslevel=2) as handle, body.open() as entries:
        handle.write("%%MatrixMarket matrix coordinate real general\n% Author values retained unchanged; GW18 subset; NOT raw counts.\n")
        handle.write("%d %d %d\n" % (nr, len(kept), kept_nz))
        for line in entries:
            handle.write(line)
    body.unlink()
    with gzip.open(out / "barcodes.tsv.gz", "wt") as handle:
        handle.write("\n".join(kept) + "\n")
    with opened(root / "features.tsv.gz") as src, gzip.open(out / "features.tsv.gz", "wt") as dst:
        feature_count = 0
        for line in src:
            dst.write(line)
            feature_count += 1
    if feature_count != nr:
        raise ValueError("Author feature count mismatch")
    with (out / "meta.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=names, delimiter="\t")
        writer.writeheader()
        writer.writerows(selected[cell] for cell in kept)
    prefixes = collections.Counter()
    with (out / "library_mapping_to_review.tsv").open("w") as handle:
        handle.write("author_prefix\trecount_library\trecount_gem\tcells\n")
        for cell in kept:
            match = re.search(r"[ACGT]{16}(?:-[1-9]\d*)?$", cell)
            if not match:
                raise ValueError("Author cell identifier needs explicit parsing: " + cell)
            prefixes[cell[:match.start()]] += 1
        for prefix, n in sorted(prefixes.items()):
            handle.write("%s\t\t\t%d\n" % (prefix, n))
    save_json(out / "subset_summary.json", dict(sample="GW18", cells=len(kept), features=nr, nonzero=kept_nz,
              expression="author values unchanged; transformation has not been inferred", prefix_counts=dict(prefixes)))
    print("GW18 author subset: %d cells; %d entries. Values unchanged." % (len(kept), kept_nz))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    sub = commands.add_parser("convert")
    sub.add_argument("output")
    sub.add_argument("--expected", type=int, default=0)
    sub.set_defaults(func=convert)
    sub = commands.add_parser("package")
    sub.add_argument("source")
    sub.add_argument("output")
    sub.set_defaults(func=package)
    sub = commands.add_parser("merge")
    sub.add_argument("manifest")
    sub.add_argument("star_root")
    sub.add_argument("output")
    sub.add_argument("--feature", choices=["Gene", "GeneFull"], required=True)
    sub.set_defaults(func=merge)
    sub = commands.add_parser("author-subset")
    sub.add_argument("directory")
    sub.add_argument("output")
    sub.set_defaults(func=author_subset)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
