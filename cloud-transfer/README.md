Windows desktop-cloud GitHub receiver
====================================

This transfers the GW18 SRA-recount software bundle through the previously
approved internal GitHub address, `172.30.100.20`. It does not transfer new
sequencing data or re-download the old 660 MB RefSeq payload.

First use
---------

Paste the contents of `install_in_git_bash.sh` into Windows Git Bash. It obtains
only the repository metadata and this small receiver script, verifies the
script SHA256, and starts it. The existing large reference repository and its
downloaded parts are left untouched.

The receiver creates two Desktop files:

* `crouch_cloud_receive.sh`: the reusable receiver.
* `crouch_cloud_receive.cmd`: double-click this after reopening the desktop VM.

Every Git call applies Git's `http.curloptResolve` setting for
`github.com:443:172.30.100.20`, HTTP/1.1, and Windows Schannel. This uses the
same authorized internal route as the earlier hosts entry, while not depending
on the Windows hosts file or DNS cache. No administrator elevation is needed.
System hosts and DNS settings, and global Git configuration, are not changed.
These command-scoped settings apply to the receiver's Git requests; they do
not configure a web browser or unrelated programs. An HTTP proxy is disabled
only for these direct GitHub requests so the address override can take effect.
HTTPS hostname and certificate verification remain enabled.

Later use
---------

Double-click `crouch_cloud_receive.cmd`, or run in Git Bash:

```bash
bash "$HOME/Desktop/crouch_cloud_receive.sh"
```

After completion, the verified archive is copied to:

```text
C:\Users\guojunping\Desktop\crouch_gw18_recount_bundle.tar.gz
```

An additional verified copy and the two cached parts are in:

```text
C:\Users\guojunping\Desktop\crouch-cloud-transfer\download\
```

If a different archive with the same name already exists on the Desktop, it
is preserved. In that case use the verified copy in the download folder.

The archive must be uploaded into the cluster `upload/Crouch/` directory,
then unpacked there, as explained in the analysis bundle's README. It is not
another RefSeq archive for `Crouch/software/`.

Retries, resume, and stop
-------------------------

* Two parts, at most 2 MiB each, fetched one Git object at a time.
* SHA256 verification per part and after reassembly, plus gzip integrity check.
* Completed verified parts are skipped on restart. An interrupted part is
  retried from the start of that small part; Git itself does not support byte
  resume of an incomplete pack.
* Up to 20 attempts per network operation, with 10/20/40/60-second backoff.
* A single attempt is stopped after 10 minutes. The receiver stops after
  30 minutes without verified part progress. Authentication, permission,
  checksum or unsupported-server errors stop instead of retrying forever.
* Stop with Ctrl+C. A `.cmd` window stays open at the end so the result is visible.
* To stop from another Git Bash window:

```bash
touch "$HOME/Desktop/crouch-cloud-transfer/STOP"
```

To resume after a STOP file:

```bash
rm -f "$HOME/Desktop/crouch-cloud-transfer/STOP"
bash "$HOME/Desktop/crouch_cloud_receive.sh"
```

`transfer.log` and `last-network.log` are kept in the Desktop
`crouch-cloud-transfer` folder. If the internal address is unavailable or the
institutional mirror has not synchronized this commit yet, the receiver will
report that problem; refreshing public DNS does not make that internal service
available. The receiver cannot make an unavailable institutional mirror work.

Release identity
----------------

Release folder: `gw18-recount-v1`

Archive: `crouch_gw18_recount_bundle.tar.gz`

Bytes: `2735965`

SHA256:
`7250b10ce5a8a1cc0f4949c5feffad9ac87a2ca2c5c309bd2e09d1134f932eaf`

The script and manifests pin this release, not whichever file happens to share
a filename in a future revision. No infinite polling or scheduled task is installed.

Local verification before publishing
------------------------------------

Tested with real local Git repositories and the actual release blobs: partial
clone filtering, one-object fetches, simulated early EOF, bounded retries,
preservation of the first completed part, resume, final SHA256, cached completion
without a new network request, STOP handling, attempt timeout and lock cleanup.
Windows Schannel, the generated `.cmd`, and the private desktop-cloud endpoint
require the user's Windows session to execute; they are not claimed as tested
from the Mac.
