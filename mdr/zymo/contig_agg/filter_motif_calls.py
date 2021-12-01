import argparse
import math

def parseArgs():
    parser=argparse.ArgumentParser(description='call meth for filtered genome positions based on given thresholds')
    parser.add_argument('-m', '--meth', type=float, required=True, help='minimum meth prob to be called methylated')
    parser.add_argument('-u', '--unmeth', type=float, required=False, help='max unmeth prob to be called unmethylated')
    parser.add_argument('-i', '--infile', type=str, required=True, help='meth info file from  megalodon_extract_barcode_methprobs')
    parser.add_argument('-o', '--outfile', type=str, required=True, help='output file')
    args=parser.parse_args()
    return args

##for testing
##infile='/mithril/Data/Nanopore/projects/methbin/zymo/contig_agg/20190809_zymo_control/20190809_zymo_control.curate_extract.csv'
##outfile='/mithril/Data/Nanopore/projects/methbin/zymo/contig_agg/20190809_zymo_control/20190809_zymo_control.curate_filter.csv'
def write_calls(called_positions, outfile):
    '''
    take called_positions list of lists
    write out csv
    '''
    with open(outfile, 'w') as f:
        for i in called_positions:
            f.write(','.join([str(x) for x in i])+'\n')

            
def filter_calls(methinfo, meth, unmeth):
    '''
    take methinfo
    return list of lists with called positions:
    [chrom, pos, strand, log meth, motif, base, call]
    NB: if unemth is not set, then anything not called as meth is unmeth
    '''
    called_positions=[]
    for i in methinfo:
        pmeth=math.exp(i[3])
        if pmeth > meth:
            i.append('m')
            called_positions.append(i)
        if unmeth is not None:
            if 1-pmeth > unmeth:
                i.append('u')
                called_positions.append(i)
        else:
            i.append('u')
            called_positions.append(i)
    return called_positions


def read_meth(infile):
    '''
    read from meth file
    returns list of lists with meth info:
    [chrom, pos, strand, log meth, motif, base]
    '''
    with open(infile) as f:
        content=f.read().split('\n')
    methinfo=[]
    for info in content[:-1]:
        i=info.split(',')
        methinfo.append([i[0], int(i[1]), i[2], float(i[3]), i[4], i[5]])
    return(methinfo)


def main(meth, unmeth, infile, outfile):
    methinfo=read_meth(infile)
    called_positions=filter_calls(methinfo, meth, unmeth)
    write_calls(called_positions, outfile)

    
if __name__ == "__main__":
    args=parseArgs()
    main(args.meth, args.unmeth, args.infile, args.outfile)
