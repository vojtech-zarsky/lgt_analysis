import re
from pathlib import Path

class Taxonomy:
    def __init__(self, sTaxdumpDir='{}/databases/taxdump'.format(str(Path.home())), sTaxdumpMod=None, bCorrections=True, bNewParse=True, bAddSupergroups=False, bCleanTaxons=False):
        self.dCorrections = {'yersinia':'629', 'morganella':'581', 'bacteria':'2'}
        self.lSupergroups = [('eukaryota','archaeplastida',['viridiplantae','rhodelphis','rhodophyta','glaucocystophyceae'])]
        self.bAddSupergroups = bAddSupergroups
        self.bCleanTaxons = bCleanTaxons
        if bNewParse:
            self.newParse(sTaxdumpDir, sTaxdumpMod)
    
    def newParse(self, sTaxdumpDir, sTaxdumpMod):
        ## read mod file
        lMods = []
        if sTaxdumpMod != None:
            for sLine in open(sTaxdumpMod):
                if sLine.strip() == '' or sLine[0] == '#':
                    continue
                lLine = sLine.strip().split('\t')
                lMods.append( lLine )
                #if lLine[0] == 'add':
                #    lAdd.append( lLine[1:] )
        self.dNamesIds = {}
        self.dIdsNames = {}
        print('reading names/ids')
        lNonscientificLines = []
        for sLine in open(sTaxdumpDir+'/names.dmp'):
            lLine = sLine.replace('\'-',' -').replace('\'','').lower().strip('	|\n').split('\t|\t') ##
            if lLine[3] == 'scientific name':
                if lLine[1] not in self.dNamesIds:
                    self.dNamesIds[lLine[1]] = []
                self.dNamesIds[lLine[1]].append(lLine[0])
            else:
                lNonscientificLines.append(lLine)
            if lLine[0] not in self.dIdsNames:
                self.dIdsNames[lLine[0]] = []
            if lLine[3] == 'scientific name':
                self.dIdsNames[lLine[0]] = [lLine[1]]+self.dIdsNames[lLine[0]]
            else:
                self.dIdsNames[lLine[0]].append(lLine[1])
        
        for lLine in lNonscientificLines:
            if lLine[1] not in self.dNamesIds:
                self.dNamesIds[lLine[1]] = [lLine[0]]
        print('reading nodes')
        self.dNodes = {}
        for sLine in open(sTaxdumpDir+'/nodes.dmp'):
            lLine = sLine.lower().strip('	|\n').split('\t|\t') ##
            self.dNodes[lLine[0]] = {'parent':lLine[1], 'rank':lLine[2]}
        
        print('reading merges')
        self.dMerges = {}
        for sLine in open(sTaxdumpDir+'/merged.dmp'):
            lLine = list(map(lambda x:x.strip(), sLine.split('|',2)[:2]))
            self.dMerges[lLine[0]] = lLine[1]
        iNewTaxId = 1
        ## add from mod
        for lMod in lMods:
            if lMod[0] == 'add':
                #self.dNamesIds = {}
                #self.dIdsNames = {}
                #self.dNodes
                lTaxonomy = []
                for sTemp in lMod[1:]:
                    sTemp = sTemp.lower().strip()
                    if sTemp.isdigit():
                        sNewTaxId = sTemp
                        sNewTaxName = 'taxid{}'.format(sTemp)
                    elif sTemp.split(' ',1)[0].isdigit():
                        sNewTaxId = sTemp.split(' ',1)[0]
                        sNewTaxName = sTemp.split(' ',1)[1]
                    elif sTemp in self.dNamesIds:
                        sNewTaxId = self.dNamesIds[sTemp][0]
                        sNewTaxName = sTemp
                    else:
                        sNewTaxId = '0{}'.format(iNewTaxId)
                        iNewTaxId += 1
                        sNewTaxName = sTemp
                    lTaxonomy.append( (sNewTaxId, sNewTaxName) )
                for iTaxLevel in range(len(lTaxonomy)):
                    (sNewTaxId, sNewTaxName) = lTaxonomy[iTaxLevel]
                    if sNewTaxId in self.dIdsNames:
                        break
                    self.dIdsNames[sNewTaxId] = [sNewTaxName]
                    self.dNamesIds[sNewTaxName] = [sNewTaxId]
                    sParent = '1'
                    if iTaxLevel+1 <  len(lTaxonomy):
                        sParent = lTaxonomy[iTaxLevel+1][0]
                    self.dNodes[sNewTaxId] = {'parent':sParent, 'rank':''}
        ## add supergroups
        if self.bAddSupergroups:
            iNewTaxId = -1
            for (sParent, sSupergroup, lChildren) in self.lSupergroups:
                ## add names/ids
                self.dIdsNames[str(iNewTaxId)] = [sSupergroup]
                print('new supergroup:', str(iNewTaxId), [sSupergroup])
                self.dNamesIds[sSupergroup] = [str(iNewTaxId)]
                self.dNodes[str(iNewTaxId)] = {'parent':self.dNamesIds[sParent][0], 'rank':''}
                for sChild in lChildren:
                    self.dNodes[self.dNamesIds[sChild][0]]['parent'] = str(iNewTaxId)
            iNewTaxId -= 1

    
    def getTaxonomy(self, sTaxName):
        lTaxonomy = []
        sTaxName = sTaxName.replace('\'-',' -').replace('\'','').lower().strip()
        if sTaxName == None or sTaxName == '':
            return lTaxonomy
        if sTaxName not in self.dNamesIds:
            sTaxName = re.sub('^uncultured','',sTaxName)
            sTaxName = re.sub('^candidatus','',sTaxName)
            sTaxName = sTaxName.strip()
            if sTaxName == "":
                return []
            if sTaxName not in self.dNamesIds:
                sTaxName = " ".join(sTaxName.split()[:2])
                if sTaxName not in self.dNamesIds:
                    sTaxName = sTaxName.split()[0]
                    if sTaxName not in self.dNamesIds:
                        return lTaxonomy
        lTaxIds = self.dNamesIds[sTaxName]
        sTaxId = lTaxIds[0]
        if len(lTaxIds) > 1 and sTaxName in self.dCorrections:
            sTaxId = self.dCorrections[sTaxName]
        
        while sTaxId != '1':
            dNode = self.dNodes[sTaxId]
            lNames = self.dIdsNames[sTaxId]
            sName = lNames[0]
            if self.bCleanTaxons:
                sName = sName.replace(';','').replace(',','')
            lTaxonomy.append((sTaxId, sName, dNode['rank']))
            sTaxId = dNode['parent']
        return lTaxonomy
    
    def getTaxonomyByTaxId(self, sTaxId):
        if sTaxId in self.dMerges:
            sTaxId = self.dMerges[sTaxId]
        lTaxonomy = []
        while sTaxId != '1':
            try:
                dNode = self.dNodes[sTaxId]
            except:
                return lTaxonomy
            lNames = self.dIdsNames[sTaxId]
            sName = lNames[0]
            if self.bCleanTaxons:
                sName = sName.replace(';','').replace(',','')
            lTaxonomy.append((sTaxId, sName, dNode['rank']))
            sTaxId = dNode['parent']
        return lTaxonomy
    
    def isTaxon(self, sTaxName):
        if sTaxName in self.dNamesIds:
            return True
        return False


