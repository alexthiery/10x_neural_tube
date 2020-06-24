#!/usr/bin/env nextflow

nextflow.preview.dsl=2


def names1 = []
new File(params.folder1).eachDir() { file->
    names1 << "$baseDir/" + params.folder1.split('/')[-1] + "/" + file.getName()   
}

def names2 = []
new File(params.folder2).eachDir() { file->
    names2 << "$baseDir/" + params.folder2.split('/')[-1] + "/" + file.getName()  
}

names = names1 + names2

def index = [
    ['THI300A1', 'cellranger_count_hh4'],
    ['THI300A3', 'cellranger_count_hh4'],
    ['THI300A4', 'cellranger_count_hh4'],
    ['THI300A6', 'cellranger_count_hh4']
    ]

def testFastq = []
for (item in index) {
    testFastq << item + names.findAll{it -> it.contains(item[0])}
}

print testFastq





