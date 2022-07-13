import sys
import Bio.KEGG.REST as REST
import io
import pandas as pd
import PyPluMA

def isMicrobe(name,MM):
    return (name in MM) and (MM[name] == 0)

def isMetabolite(name,MM):
    return (name in MM) and (MM[name] == 1)

def to_df(result):
     return pd.read_table(io.StringIO(result), header=None)

def findPath(microbe, metabolite, realmicrobe, realmetabolite):
    microbe = microbe[1:len(microbe)-1]
    metabolite = metabolite[1:len(metabolite)-1]
    result = REST.kegg_list("pathway", microbe).read()
    result_frame = to_df(result)
    #print(microbe)
    for i in range(len(result_frame[0])):
        path = result_frame[0][i]
        #print("PATH: "+path)
        path = path[path.find('path:'+microbe)+len(microbe)+5:]
        path = "map"+path
        comps = REST.kegg_link("compound", path).read()
        #print("COMPS: "+comps)
        #print("LENGTH: "+str(len(comps)))
        try:
         allcomps = to_df(comps)
         for j in range(len(allcomps[0])):
            mycomp = allcomps[1][j]
            mycomp = mycomp[mycomp.find('cpd:')+4:]
            if (metabolite == mycomp):
                print(realmicrobe+","+realmetabolite+","+path+":"+result_frame[1][i])
            #else:
            #    print(metabolite+","+mycomp)
        except Exception:
             continue

def findMetabolicPath(metabolite1, metabolite2, realmetabolite1, realmetabolite2):
    metabolite1 = metabolite1[1:len(metabolite1)-1]
    metabolite2 = metabolite2[1:len(metabolite2)-1]
    #print(metabolite1)
    result = REST.kegg_link("pathway", metabolite1).read()
    try:
     result_frame = to_df(result)
     #print(microbe)
     for i in range(len(result_frame[0])):
        path = result_frame[1][i]
        path = path[path.find('path:')+5:]
        #print("PATH: "+path)
        #path = "map"+path
        comps = REST.kegg_link("compound", path).read()
        #print("COMPS: "+comps)
        #print("LENGTH: "+str(len(comps)))
        try:
         allcomps = to_df(comps)
         for j in range(len(allcomps[0])):
            mycomp = allcomps[1][j]
            mycomp = mycomp[mycomp.find('cpd:')+4:]
            if (metabolite2 == mycomp):
                pathinfo = to_df(REST.kegg_find("pathway", path).read())
                print(realmetabolite1+","+realmetabolite2+","+path+":"+pathinfo[1][0])
            #else:
            #    print(metabolite+","+mycomp)
        except Exception:
             print("TO_DF NOTHING")
             continue
    except Exception:
        print("TO_DF NOTHING")
        return


class KEGGMultiPlugin:
    def input(self, inputfile):
        params = dict()
        infile = open(inputfile, 'r')
        for line in infile:
            contents = line.strip().split('\t')
            params[contents[0]] = contents[1]
        self.csvfile = open(PyPluMA.prefix()+"/"+params['csvfile'], 'r')
        self.keggfile = open(PyPluMA.prefix()+"/"+params['keggfile'], 'r')

    def run(self):
        microbes = dict()
        metabolites=dict()
        MM=dict()
        for line in self.keggfile:
            contents =line.strip().split(',')
            MM[contents[0]] = int(contents[2])
            if (MM[contents[0]] == 0):
                microbes[contents[0]] = contents[1]
            else:
                metabolites[contents[0]] = contents[1]

        header = self.csvfile.readline().strip()
        headercontents = header.split(',')
        
        posx = 1
        for line in self.csvfile:
            contents =line.strip().split(',')
            x = contents[0]
            pos = posx
            while (pos < len(headercontents)):
                y = headercontents[pos]
                value = float(contents[pos])

                if (value > 0):
                   if (isMicrobe(x,MM) and isMetabolite(y,MM)):
                      print("Exploring: "+x+" and "+y)
                      findPath(microbes[x], metabolites[y], x, y)
                   elif (isMicrobe(y,MM) and isMetabolite(x,MM)):
                      print("Exploring: "+y+" and "+x)
                      findPath(microbes[y], metabolites[x], y, x)
                   elif (isMetabolite(x,MM) and isMetabolite(y,MM) and x != y):
                      print("Exploring: "+x+" and "+y)
                      findMetabolicPath(metabolites[x], metabolites[y], x, y)

                pos += 1
            posx += 1
    def output(self, outfile):
       pass
