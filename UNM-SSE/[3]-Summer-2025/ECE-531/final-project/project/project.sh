#!/usr/bin/env bash
set -euo pipefail

# Always run relative to where this script lives
DIR="$(cd "$(dirname "$0")" && pwd)"

# --- Start/verify server.py on host ---
if ! pgrep -f "[p]ython3 $DIR/server.py" >/dev/null ; then
  echo "[HOST] Starting server.py..."
  nohup python3 "$DIR/server.py" > "$DIR/server.log" 2>&1 &
  # Wait until port 8000 answers
  echo -n "[HOST] Waiting for server on :8000"
  for i in {1..30}; do
    if curl -fsS http://127.0.0.1:8000/ >/dev/null 2>&1 || \
       curl -fsS http://127.0.0.1:8000/status >/dev/null 2>&1; then
      echo "  OK"
      break
    fi
    echo -n "."
    sleep 0.5
  done
else
  echo "[HOST] server.py already running."
fi

# --- Launch QEMU (all paths from $DIR) ---
echo "[HOST] Starting QEMU (serial console)..."
exec qemu-system-arm \
  -M versatilepb \
  -kernel "$DIR/zImage" \
  -dtb "$DIR/versatile-pb.dtb" \
  -drive file="$DIR/rootfs.ext2",if=scsi,format=raw \
  -append "root=/dev/sda console=ttyAMA0,115200" \
  -nographic \
  -audio none \
  -net nic,model=rtl8139 \
  -net user,hostfwd=tcp::2222-:22
