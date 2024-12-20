import sys
import os
import re
import shutil
import argparse
import logging
from pybedtools import BedTool
from qualiffy import phasing

logger = logging.getLogger(__name__)

def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        usage="%(prog)s [OPTION] [FILE]...",
        description="Print phase block stats from a bed file of haplotype marker locations"
    )
    parser.add_argument(
        "-v", "--version", action="version",
        version = f"{parser.prog} version 0.1.0"
    )
    parser.add_argument('-b', '--bedfile', required=True, default=None, help='bed file with locations of haplotype markers and fourth column haplotype names')
    parser.add_argument('--shortnum', type=int, required=False, default=1, help='number of switches allowed in a region before a short block becomes a long block')
    parser.add_argument('--shortlimit', type=int, required=False, default=200, help='number of bases stretch before a short switch becomes a long block')
    parser.add_argument('-p', '--prefix', type=str, required=False, default='out', help='prefix for output phase block file and log file')
    parser.add_argument('--breakatgaps', action='store_true', required=False, default=False, help='break phase blocks at gaps in bed file')
    parser.add_argument('--debug', action='store_true', required=False, default=False, help='print verbose output')

    return parser

def parse_arguments(args):
    parser = init_argparse()
    args = parser.parse_args(args)

    return args

def main() -> None:

    args = parse_arguments(sys.argv[1:])

    logfile = args.prefix + ".log"
    logformat = '%(asctime)s %(message)s'
    if args.debug:
        logging.basicConfig(filename=logfile, level=logging.DEBUG, format=logformat)
        logger.info('Logging verbose output for debugging.')
    else:
        logging.basicConfig(filename=logfile, level=logging.INFO, format=logformat)

    logger.info("Calling find_phase_blocks_from_marker_bed")
    phasing.find_phase_blocks_from_marker_bed(args.bedfile, args.shortnum, args.shortlimit, args.breakatgaps)

if __name__ == "__main__":
    main()
