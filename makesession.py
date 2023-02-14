#!/Library/Frameworks/Python.framework/Versions/3.9/bin/python3

##making igv sessions is too long
#colours - rgb strings
outfile="test.xml"

#file path prefix
prefix = """https://www.bioinfo.ieo.eu/GNgroup/IEO/gmandana/project/drb3/"""
#num tracks
annotation = "https://www.bioinfo.ieo.eu/GNgroup/IEO/gmandana/annotation/gencodev37_annotation.sorted.gtf"
#xml header
xml_header = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="hg38" hasGeneTrack="false" hasSequenceTrack="true" locus="chr13:79357287-79455472" nextAutoscaleGroup="16" version="8">"""

#colour vector - same length as bw vector
#col_vector = []

juncs_filenames = []
juncs_names = []
juncs_url = [prefix + i for i in juncs_filenames]

bw_filenames = []
track_names = []
file_url = [prefix + i for i in bw_filenames]

##write header
with open(outfile, 'w') as f:
    f.write(xml_header)

##write resources

with open(outfile, 'a') as f:
    f.write("\t<Resources>\n")
    for file in file_url:
        f.write("\t\t<Resource path=\"{0}\"/>\n".format(file))
    for file in juncs_url:
        f.write("\t\t<Resource path=\"{0}\"/>\n".format(file))
    f.write("\t\t<Resource path=\"{0:s}\"/>\n".format(annotation))
    f.write("\t</Resources>\n")


##tracks

with open(outfile, 'a') as f:
    f.write("<Panel height=\"1920\" name=\"DataPanel\" width=\"2112\">\n")
    for i in range(len(bw_filenames)):
        f.write("\t\t<Track attributeKey=\"{0:s}\" autoScale=\"false\" autoscaleGroup=\"11\" clazz=\"org.broad.igv.track.DataSourceTrack\" fontSize=\"10\" id=\"{1:s}\" name=\"{2:s}\" renderer=\"BAR_CHART\" visible=\"true\" windowFunction=\"mean\">\n".format(bw_filenames[i], file_url[i], track_names[i]))
        f.write("\t\t\t<DataRange baseline=\"0.0\" drawBaseline=\"true\" flipAxis=\"false\" maximum=\"53.7714\" minimum=\"0.0\" type=\"LINEAR\"/>\n")
        f.write("\t\t</Track>\n")
    f.write("</Panel>\n")

##feature panel

with open(outfile, 'a') as f:
    f.write("<Panel height=\"1920\" name=\"FeaturePanel\" width=\"2112\">")
    f.write("\t\t<Track attributeKey=\"Reference sequence\" clazz=\"org.broad.igv.track.SequenceTrack\" fontSize=\"10\" id=\"Reference sequence\" name=\"Reference sequence\" sequenceTranslationStrandValue=\"POSITIVE\" shouldShowTranslation=\"false\" visible=\"true\"/>")
    f.write("\t\t<Track attributeKey=\"gencodev37_annotation.sorted.gtf\" clazz=\"org.broad.igv.track.FeatureTrack\" displayMode=\"EXPANDED\" featureVisibilityWindow=\"100000\" fontSize=\"10\" groupByStrand=\"false\" id=\"https://www.bioinfo.ieo.eu/GNgroup/IEO/gmandana/annotation/gencodev37_annotation.sorted.gtf\" name=\"gencodev37_annotation.sorted.gtf\" visible=\"true\"/>")
    for j in range(len(juncs_names)):
        f.write("\t\t<Track attributeKey=\"{0:s}}\" clazz=\"org.broad.igv.track.FeatureTrack\" colorScale=\"ContinuousColorScale;0.0;1106.0;255,255,255;0,0,178\" fontSize=\"10\" groupByStrand=\"false\" height=\"60\" id=\"{1:s}\" name=\"{2:s}\" visible=\"true\"/>".format(juncs_filenames[j], juncs_url[j], juncs_names[j]))
    f.write("</Panel>")


##footer stuff
with open(outfile, 'a') as f:
    f.write("\t<PanelLayout dividerFractions=\"0.5932377049180327\"/>\n")
    f.write("\t<HiddenAttributes>\n")
    f.write("\t\t<Attribute name=\"DATA FILE\"/>\n")
    f.write("\t\t<Attribute name=\"DATA TYPE\"/>\n")
    f.write("\t\t<Attribute name=\"NAME\"/>\n")
    f.write("\t</HiddenAttributes>\n")

##close

with open(outfile, 'a') as f:
    f.write("</Session>")





