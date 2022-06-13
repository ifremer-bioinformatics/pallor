#!/usr/bin/env python3

import argparse
import os


def getArgs():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-f', dest="folder", type=str, required=True, help='Single copy BUSCO sequences folder')
    parser.add_argument('-d', dest="dbodb", type=str, required=True, help='BUSCO ODB_10 database name')
    parser.add_argument('-s', dest="specie", type=str, required=True, help='Specie ID')

    arg = parser.parse_args()

    return arg


def main(args):
    cat_faa_sg = open(args.specie + '.sg.faa', 'w')
    busco_single_copy_folder = os.path.join(args.folder, 'run_' + args.dbodb, 'busco_sequences/single_copy_busco_sequences')
    for file in os.listdir(busco_single_copy_folder):
        if file.endswith('.faa'):
            faa = open(os.path.join(busco_single_copy_folder, file), 'r')
            odb_id = os.path.splitext(file)[0]
            for line in faa:
                if line.startswith('>'):
                    ctg_start_end = line.split()[0].replace('>', '')
                    cat_faa_sg.write('>' + odb_id + ':' + args.specie + ':' + ctg_start_end + '\n')
                else:
                    cat_faa_sg.write(line)


if __name__ == '__main__':
    args = getArgs()
    main(args)
