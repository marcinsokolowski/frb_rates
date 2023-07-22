#!/bin/bash

options=""
if [[ -n "$1" && "$1" != "-" ]]; then
   options="$1"
fi

echo "python ~/github/frb_rates/python/plot_frb_rates.py --grid --legend_with_curves ${options}"
python ~/github/frb_rates/python/plot_frb_rates.py --grid --legend_with_curves ${options}
