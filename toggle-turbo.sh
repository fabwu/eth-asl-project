#!/bin/bash

TURBOFILE="/sys/devices/system/cpu/intel_pstate/no_turbo"

turbo_disabled=$(cat "$TURBOFILE")

if [ "$turbo_disabled" = "1" ];
then
  echo '0' | sudo tee -a "$TURBOFILE" > /dev/null
  echo "Turbo Enabled"
else
  echo '1' | sudo tee -a "$TURBOFILE" > /dev/null
  echo "Turbo Disabled"
fi
