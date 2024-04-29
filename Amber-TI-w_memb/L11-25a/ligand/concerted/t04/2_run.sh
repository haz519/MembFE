#!/usr/bin/sh

# arg1 is the step

if [ $# -lt 1 ]; then
  echo "Usage: $0 step (e.g., 2_run.sh 1)"
  exit 1
fi

case "$1" in
  [0-9]*)
    start=$1
  ;;
  *)
    echo "start must be a number"
    exit 1;
  ;;
esac

for w in $(cat mbar.dat); do
    cd $w
    ../submit.sh $1
    cd ..
done
