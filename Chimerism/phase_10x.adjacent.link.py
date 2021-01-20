#usr/bin/env python3
import sys
import copy

if len(sys.argv) < 4:
    sys.exit("python3 {} *.snp.reads win step".format(sys.argv[0]))

def snps2allele(snps):
    allele = ""
    for a,b in snps:
        allele += b
    return allele

class TreeNode:
    def __init__(self, val):
        self.val = val
        self.left, self.right = None, None


def dfs(node, result, tmp):
    if node is None:
        return

    tmp.append(node)
    # 这里需要用拷贝而不是用 = 赋值，也可以遍历赋值
    tmp1 = copy.deepcopy(tmp)

    if node.left is None and node.right is None:
        result.append([i.val for i in tmp])
        return

    if node.left is not None:
        dfs(node.left, result, tmp)
    # 遍历右子树需要带上不同的变量，否则左子树的tmp和右子树的tmp都指向一块内存
    if node.right is not None:
        dfs(node.right, result, tmp1)

def is_odd(dd):
    if dd % 2 == 0:
        return True
    else:
        return False

def depth_adjacent(site1, site2, records):
    (pos1, base1) = site1
    (pos2, base2) = site2
    # calculate how many 10x read support
    barcodes = {}
    # reads mapped to pos1
    for l in records[pos1]:
        tmp = l.split("\t")
        read_id = tmp[-1]
        bar = read_id.split("_")[-1]
        if tmp[4] == "ref":
            tmp_site = (pos1, tmp[2])
        else:
            tmp_site = (pos1, tmp[3])
        if bar not in barcodes.keys():
            barcodes[bar] = [tmp_site,]
        else:
            if tmp_site not in barcodes[bar]:
                barcodes[bar].append(tmp_site)
    # reads mapped to pos2
    for l in records[pos2]:
        tmp = l.split("\t")
        read_id = tmp[-1]
        bar = read_id.split("_")[-1]
        if tmp[4] == "ref":
            tmp_site = (pos2, tmp[2])
        else:
            tmp_site = (pos2, tmp[3])
        if bar not in barcodes.keys():
            barcodes[bar] = [tmp_site,]
        else:
            if tmp_site not in barcodes[bar]:
                barcodes[bar].append(tmp_site)


    read_support = 0
    for b in barcodes.keys():
        if site1 in barcodes[b] and site2 in barcodes[b]:
            read_support += 1
    return read_support


def makeup_allele(allele, pre, s):
    (pos1, base1) = pre
    (pos2, base2) = s
    if len(allele) == 0:
        return base1 + base2
    else:
        return allele + base2



def main():
    molecular = {}
    snps = set()
    records = {}
    ref_info = []
    alt_info = []
    # store all records in dict, key=position; val=list of reads
    with open(sys.argv[1], 'r') as fh:
        for line in fh:
            line = line.strip()
            tmp = line.split("\t")
            pos = int(tmp[1])
            snps.add(pos)
            if (pos, tmp[2]) not in ref_info:
                ref_info.append((pos, tmp[2]))
            if (pos, tmp[3]) not in alt_info:
                alt_info.append((pos, tmp[3]))
            if pos not in records.keys():
                records[pos] = [line,]
            else:
                records[pos].append(line)

    # deal by windows and slide by step
    poses = sorted(list(snps))
    win = int(sys.argv[2])
    step = int(sys.argv[3])
    start = 0
    end = start + win
    while end <= len(poses):
        pos_info = ",".join([str(i) for i in poses[start:end]])
        tmp_ref_info = ref_info[start:end]
        tmp_alt_info = alt_info[start:end]
        refseq = snps2allele(tmp_ref_info)
        #list1 = ['a', 'b', 'c', 'd']
        #list2 = ['e', 'f', 'g', 'h']

        nodes = {}
        # store the node information
        level = len(tmp_ref_info) + 1
        for i in range(level):
            if i == 0:
                nodes[0] = TreeNode("0")
            else:
                this_level_min = 2 ** i - 1
                this_level_max = 2 ** (i + 1) - 2
                #print("{}\t{}".format(this_level_min, this_level_max))
                for j in range(this_level_min, this_level_max+1):
                    if is_odd(j):
                        nodes[j] = TreeNode(tmp_alt_info[i-1])
                    else:
                        nodes[j] = TreeNode(tmp_ref_info[i-1])

        # store node relationship
        """
         artifical           node1("0")
                        /                \
         allele1    node3(ref)           node4(alt)
                    /      \             /         \
         allele2 node7(ref) node8(alt) node9(ref)  node10(alt)
        """
        for i in range(level-1):
            this_level_min = 2 ** i - 1
            this_level_max = 2 ** (i + 1) - 2
            for j in range(this_level_min, this_level_max+1):
                nodes[j].left = nodes[2*j+1]
                nodes[j].right = nodes[2*j+2]


        paths = []
        tmp = [] # be carefully
        results = {}
        results_subdepth = {}
        # get all potential paths, traveling binary tree graph
        dfs(nodes[0], paths, tmp)
        for path in paths:
            path.remove('0') # remove the first element "0"
            allele = ""
            allele_depth = 0
            sub_depth = []
            for s in path:
                if len(allele) == 0:
                    allele = s[1]   # s = (pos, base)
                    pre = s
                else:
                    this_depth = depth_adjacent(pre, s, records)
                    if this_depth > 0:
                        allele = makeup_allele(allele, pre, s)
                        allele_depth += this_depth
                        sub_depth.append(this_depth)
                        pre = s
                    else:
                        # if this two-adjacent-snp not support by any reads, remove this
                        # path
                        break
                if len(allele) == win:
                    results[allele] = allele_depth
                    results_subdepth[allele] = sub_depth
        # sort alleles by depth
        sorted_alleles = sorted(results, key=lambda k: (results[k], k), reverse=True)
        if len(sorted_alleles) > 0:
            print(">" + pos_info + "\t" + refseq + "\t" + str(len(sorted_alleles)))
        for k in sorted_alleles:
            print(k + "\t" + str(results[k]) + "\t" + ",".join([ str(i)  for i in results_subdepth[k] ]))

        start += step
        end = start + win


if __name__ == '__main__':
    main()


