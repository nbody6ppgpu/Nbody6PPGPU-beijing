#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 

import sys
import time
import re

def main():
    timestamp = None
    for line in sys.stdin.readlines():
        line = line.rstrip('\n')
        # line example: ./3node-12mpi-12gpu-OMPclose-visible-OMPPlace/N1m.10472543.out-   0 12 1000000    797.95063     52.31    730.16    172.58     29.70     65.19     22.36     41.91      0.00     99.04     25.11     70.83     25.82     20.16     14.62     15.24      0.02    238.60     -0.00    120.16     71.87     11.26    150.02      0.02      0.00      0.00      0.00      0.00      9.51      0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00     10.73      3.81      0.02      0.00      0.00      0.00      0.00  3.59601E+11  8.85213E+11  0.00E+00  0.00E+00         0         0
        # 将line用空格分割
        parts = line.split()
        # parts[0] example: ./3node-12mpi-12gpu-OMPclose-visible-OMPPlace/N1m.10472543.out-
        # parts[0]分割出来几个值: jobID (10472543), #node (3), #mpi (12), #GPU (12), 其他参数
        jobID = parts[0].split('/')[-1].split('.')[1]
        # 使用re搜索node,mpi,gpu
        try:
            node = re.search(r'(\d+)node', parts[0]).group(1)
        except AttributeError:
            node = '?'
        try:
            mpi = re.search(r'(\d+)mpi', parts[0]).group(1)
        except AttributeError:
            mpi = '?'
        try:
            gpu = re.search(r'(\d+)gpu', parts[0]).group(1)
        except AttributeError:
            gpu = '?'
        other_labels = parts[0].split('gpu')[-1].split('/')[0].split('-')[1:]

        # 整合数据：先是parts1，然后jobid, node, mpi, gpu, others，使用逗号分隔，输出为文本
        sys.stdout.write(','.join(parts[1:] + [jobID, node, mpi, gpu] ) + ',' + '-'.join(other_labels) + '\n')

if __name__ == '__main__':
    main()
