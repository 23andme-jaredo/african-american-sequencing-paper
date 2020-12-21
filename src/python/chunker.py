#! /usr/bin/env python

import argparse
import logging
import hashlib
import subprocess

import seqUtils


def main(args):
    tasks = []
    logging.info("Processing %s"%args.vcf)
    for c,a,b in seqUtils.make_chunks(args.w,args.b,args.vcf):
        region = "%s:%d-%d"%(c,a,b)
        print(region)

if __name__ == "__main__":    
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf', type=str)
    parser.add_argument('-w', type=int,default=1000000)
    parser.add_argument('-b', type=int,default=200000)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    main(args)
