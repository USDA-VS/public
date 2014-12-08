#! /usr/bin/env python

# Finalized by Kaity Brien on 2014/11/06
# Python script to format kraken report and get hierarchical list of TaxIDs

from sys import argv
import re
import os
from operator import eq

root= os.getcwd()

# Format kraken report to get rid of spaces
# open file
script, filename, organism = argv
print("script is " + script)
print("filename is " + filename)
print("organism is " + organism)
file= open(filename)
#re.sub(patternToSearch, replacementDesired, stringToSearch, optionalCount, optionalFlags)
fileNameToSave= re.sub(".txt", "", filename)

# open tempfile to save formatted content to
formattedFileName= fileNameToSave + "FormattedFile.txt"
formattedFile= open(formattedFileName, "w+")


#Format output of kraken-report to prepare for processing
for line in file:
    #get rid of spaces at beginning of file
    replaceFirstSpace= re.sub(" ", "", line, count=1)
    #get rid of two spaces separating taxonomy categories and replace with tab
    replaceSpaceWTabs= re.sub("\s\s", "\t", replaceFirstSpace)
    #get rid of last spaces
    #replaceLastSpaces= re.sub(" ", "", replaceSpaceWTabs)
    replaceLastSpaces=re.sub(" ", "_", replaceSpaceWTabs)
    #replaceLeadingUnderscore=re.sub("\t_", "", replaceLastSpaces)
    formattedFile.write(replaceLastSpaces)
formattedFile.close()
file.close()

#save each line of file as an arrayList
fileToParse= open(formattedFileName)

# Initialize matrix for storing hierarchy
hierarchyMatrix= []
# Add each line in file as a list to the list of matrices (list with each line in file as a list)
for line in fileToParse:
    newLine= line.rstrip('\n').split('\t')
    hierarchyMatrix.append(newLine)

# Fill in holes (empty strings) in matrix
previousRow= []
lengthLongestRow = 0
for list in hierarchyMatrix:
    currentRow= list
    lengthOfCurrentRow= len(currentRow)
    lengthOfPreviousRow= len(previousRow)
    if(lengthOfCurrentRow > lengthLongestRow):
        lengthLongestRow = lengthOfCurrentRow
    counter=0
    while((counter < lengthOfCurrentRow) and (counter< lengthOfPreviousRow)):
        if (previousRow != []):
            if(previousRow[counter] is not None):  #checks for nullness
                if(currentRow[counter]==''):
                    currentRow[counter]=previousRow[counter] 
        counter+=1
    previousRow = currentRow

# Make folder for files that will be generated
currentPath= os.getcwd()
folderName= currentPath + "/hierarchicalClustering"
os.mkdir(folderName)
os.chdir(folderName)
currentPath=os.getcwd()

# Cluster based on different levels of hierarchy
# Start of column 6 of 1st list

columnCounter= 5
rowCounter= 0
thisRow= hierarchyMatrix[rowCounter]
listOfNames=[]
# Haven't exhausted columns
while (columnCounter < lengthLongestRow):
    # Haven't exhausted the rows
    columnAsString= str(columnCounter)
    newDirName= currentPath + "/column" + columnAsString
    newPath= currentPath + "/column" + columnAsString
    os.mkdir(newDirName)
    os.chdir(newPath)
    if (rowCounter==len(hierarchyMatrix)):
        rowCounter=0
    while (rowCounter < len(hierarchyMatrix)):
        thisRow=hierarchyMatrix[rowCounter]
        lengthCurrentRow = len(thisRow)
        # Check if aren't going to reference beyond list index
        if(columnCounter < lengthCurrentRow):
            name= thisRow[columnCounter]
            formattedName= re.sub("/", "-", name)
            # Can change this conditional to require greater number of reads to be included in taxonIDs to search
            if (int(thisRow[2]) > 0):  # Number of reads assigned directly to this taxon >0
                if (formattedName not in listOfNames):
                    listOfNames.append(formattedName)
                    nameOfOutfile= formattedName + "ListOfTaxIDs"
                    outfile= open(nameOfOutfile, "w")
                    taxidToWrite= hierarchyMatrix[rowCounter][4]
                    outfile.write(taxidToWrite + '\n')
                elif (formattedName in listOfNames):
                    taxidToWrite= hierarchyMatrix[rowCounter][4]
                    outfile.write(taxidToWrite + '\n')
        rowCounter+=1
    outfile.close()
    columnCounter+=1
#print("Reached end of this script")                







    

