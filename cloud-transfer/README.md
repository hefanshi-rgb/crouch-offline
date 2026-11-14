# Metadata relay

The completed large software/reference/recount transfers were retired.

Run `python3 tools/prepare_author_meta_on_mac.py` on the Mac after merging the cleanup commit. Select the author complete meta.tsv. The helper prepares `author-meta-v1/`, a hash-pinned `crouch_author_meta_receive.sh` here, and a Windows Git Bash bootstrap text on the Mac Desktop. Commit the generated new files and push using GitHub Desktop.

The Windows receiver uses the approved command-scoped GitHub route, verifies every 2 MiB part and the assembled archive, and retains verified progress on failure. Ctrl+C or its STOP file halts the transfer. It does not edit hosts, DNS or global Git configuration.
