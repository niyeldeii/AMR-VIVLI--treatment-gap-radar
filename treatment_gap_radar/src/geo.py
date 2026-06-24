"""ISO-3 -> continent and World Bank income tier (approximate, for modelling features).

Income tiers are coarse World Bank groupings (H=high, UM=upper-middle, LM=lower-middle, L=low);
used only as a model feature and clearly labelled approximate.
"""
CONTINENT = {
 "ARE":"Asia","ARG":"South America","AUS":"Oceania","AUT":"Europe","BEL":"Europe","BGR":"Europe",
 "BLR":"Europe","BRA":"South America","CAN":"North America","CHE":"Europe","CHL":"South America",
 "CHN":"Asia","CIV":"Africa","CMR":"Africa","COL":"South America","CRI":"North America","CZE":"Europe",
 "DEU":"Europe","DNK":"Europe","DOM":"North America","ECU":"South America","EGY":"Africa","ESP":"Europe",
 "EST":"Europe","FIN":"Europe","FRA":"Europe","GBR":"Europe","GHA":"Africa","GRC":"Europe",
 "GTM":"North America","HKG":"Asia","HND":"North America","HRV":"Europe","HUN":"Europe","IDN":"Asia",
 "IND":"Asia","IRL":"Europe","ISR":"Asia","ITA":"Europe","JAM":"North America","JOR":"Asia","JPN":"Asia",
 "KEN":"Africa","KHM":"Asia","KOR":"Asia","KWT":"Asia","LBN":"Asia","LTU":"Europe","LVA":"Europe",
 "MAR":"Africa","MEX":"North America","MUS":"Africa","MWI":"Africa","MYS":"Asia","NAM":"Africa",
 "NGA":"Africa","NIC":"North America","NLD":"Europe","NOR":"Europe","NZL":"Oceania","OMN":"Asia",
 "PAK":"Asia","PAN":"North America","PHL":"Asia","POL":"Europe","PRI":"North America","PRT":"Europe",
 "QAT":"Asia","ROU":"Europe","RUS":"Europe","SAU":"Asia","SGP":"Asia","SLV":"North America",
 "SRB":"Europe","SVK":"Europe","SVN":"Europe","SWE":"Europe","THA":"Asia","TUN":"Africa","TWN":"Asia",
 "UGA":"Africa","UKR":"Europe","USA":"North America","VEN":"South America","VNM":"Asia","ZAF":"Africa",
}
INCOME = {
 "ARE":"H","ARG":"UM","AUS":"H","AUT":"H","BEL":"H","BGR":"UM","BLR":"UM","BRA":"UM","CAN":"H","CHE":"H",
 "CHL":"H","CHN":"UM","CIV":"LM","CMR":"LM","COL":"UM","CRI":"UM","CZE":"H","DEU":"H","DNK":"H","DOM":"UM",
 "ECU":"UM","EGY":"LM","ESP":"H","EST":"H","FIN":"H","FRA":"H","GBR":"H","GHA":"LM","GRC":"H","GTM":"UM",
 "HKG":"H","HND":"LM","HRV":"H","HUN":"H","IDN":"UM","IND":"LM","IRL":"H","ISR":"H","ITA":"H","JAM":"UM",
 "JOR":"LM","JPN":"H","KEN":"LM","KHM":"LM","KOR":"H","KWT":"H","LBN":"LM","LTU":"H","LVA":"H","MAR":"LM",
 "MEX":"UM","MUS":"UM","MWI":"L","MYS":"UM","NAM":"UM","NGA":"LM","NIC":"LM","NLD":"H","NOR":"H","NZL":"H",
 "OMN":"H","PAK":"LM","PAN":"H","PHL":"LM","POL":"H","PRI":"H","PRT":"H","QAT":"H","ROU":"H","RUS":"UM",
 "SAU":"H","SGP":"H","SLV":"LM","SRB":"UM","SVK":"H","SVN":"H","SWE":"H","THA":"UM","TUN":"LM","TWN":"H",
 "UGA":"L","UKR":"LM","USA":"H","VEN":"UM","VNM":"LM","ZAF":"UM",
}
