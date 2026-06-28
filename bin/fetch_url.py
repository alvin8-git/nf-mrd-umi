#!/usr/bin/env python3
"""Resumable HTTP downloader (curl/wget are blocked in this env).

Resumes a partial file via a Range request and retries on transient errors, so a
27 GB download over throttled/flaky egress survives drops instead of restarting.
Verifies the final size against the server's reported total.

Usage: fetch_url.py URL DEST [URL DEST ...]
"""
import os, sys, time, urllib.request, urllib.error

def total_size(url):
    req = urllib.request.Request(url, headers={"Range": "bytes=0-0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        cr = r.headers.get("Content-Range")
        return int(cr.split("/")[-1]) if cr else None

def fetch(url, dest, retries=100):
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    total = total_size(url)
    for attempt in range(retries):
        have = os.path.getsize(dest) if os.path.exists(dest) else 0
        if total and have >= total:
            print(f"[done] {dest} ({have/1e9:.2f} GB)"); return
        hdrs = {"Range": f"bytes={have}-"} if have else {}
        try:
            req = urllib.request.Request(url, headers=hdrs)
            with urllib.request.urlopen(req, timeout=120) as r, open(dest, "ab") as f:
                t0, last = time.time(), have
                while True:
                    chunk = r.read(1 << 20)  # 1 MB
                    if not chunk:
                        break
                    f.write(chunk); have += len(chunk)
                    if time.time() - t0 > 30:
                        rate = (have - last) / (time.time() - t0) / 1e6
                        pct = f"{100*have/total:.1f}%" if total else "?"
                        print(f"[{os.path.basename(dest)}] {have/1e9:.2f} GB {pct} ~{rate:.1f} MB/s", flush=True)
                        t0, last = time.time(), have
        except (urllib.error.URLError, ConnectionError, TimeoutError) as e:
            wait = min(60, 2 ** min(attempt, 6))
            print(f"[retry {attempt+1}] {os.path.basename(dest)} @ {have/1e9:.2f} GB: {e} (sleep {wait}s)", flush=True)
            time.sleep(wait); continue
        # stream ended; loop re-checks completeness (handles silent truncation)
    final = os.path.getsize(dest)
    if total and final < total:
        sys.exit(f"INCOMPLETE {dest}: {final}/{total} after {retries} attempts")
    print(f"[done] {dest} ({final/1e9:.2f} GB)")

if __name__ == "__main__":
    a = sys.argv[1:]
    if len(a) < 2 or len(a) % 2:
        sys.exit("usage: fetch_url.py URL DEST [URL DEST ...]")
    for i in range(0, len(a), 2):
        print(f"=== fetching {a[i]} -> {a[i+1]} ===", flush=True)
        fetch(a[i], a[i+1])
    print("ALL DOWNLOADS COMPLETE")
