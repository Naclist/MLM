# -*- coding: utf-8 -*-
"""
MLM Pipeline - Multi-Level-Migration Analysis
Author: Naclist
Version: 0.1
Date: 2024/11/29

This script is a submodule of Multi-Level-Migration (MLM).
"""
import sys


def parse_nex(nexus, blank_value='Unknown'):
    with open(nexus, 'r') as nxin:
        with open(nexus+'.meta', 'w') as nout:
            for i in nxin:
                if i.strip().startswith('Tree'):
                    for j in i.strip().split(']'):
                        part = j.split('(')[-1].replace('"', '').replace(')', '').replace(',', '').split('[&state=')
                        if len(part) > 1:
                            print(f"{part[0].split(':')[0]}\t{part[1]}", file=nout)
                        else:
                            print(f"{part[0].split(':')[0]}\t{blank_value}", file=nout)


if len(sys.argv) < 2:
    print('USAGE: [parse_annotated_nex] [blank value (Default: Unknown)]')
else:
    nexus = sys.argv[1]
    if len(sys.argv) > 2:
        parse_nex(nexus, sys.argv[2])
    else:
        parse_nex(nexus, )