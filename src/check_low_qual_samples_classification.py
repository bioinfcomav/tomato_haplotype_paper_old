import pandas

classifications = pandas.read_csv('/home/david/Disco/ultimate_paper/classifications.csv', index_col='sample')
with open('/home/david/Disco/ultimate_paper/source_data/low_quality_samples_085.txt') as low_qual_samples_file:
    low_qual_samples = [line.rstrip() for line in low_qual_samples_file]

for sample in low_qual_samples:
    print(sample.lower(),  end='\t')
    try:
        classification = classifications['rank2'][sample.lower()]
    except KeyError:
        classification = 'None'
    print(classification)