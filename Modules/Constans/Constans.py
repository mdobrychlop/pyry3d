#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# 
# www.genesilico.pl 
#

#creates ranked 3D models of macromoleular complexes 
#based on experimental restraints and a whole complex shape.


__author__ = "Joanna M. Kasprzak"
__copyright__ = "Copyright 2010, The PyRy3D Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Joanna Kasprzak"
__email__ = "jkasp@amu.edu.pl"
__status__ = "Prototype"

#######################  BioPython version   ############################
        
#######################   path settings  ############################


#######################   Aminoacids and nucleotides   ###################
#contains both amino acids and nucleotides
RESNAMES = {"ALA": "A", "ARG": "R", "ASP": "D", "ASN": "N", "CYS": "C",\
              "GLU": "E", "GLY": "G", "GLN": "Q", "HIS": "H", \
              "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "MSE": "M",\
              "PHE": "F", "PRO": "P", "SER": "S", "THR": "T",\
              "TRP": "W", "TYR": "Y", "VAL": "V",
              "CYT": "C", "THY": "T",\
              "GUA": "G","ADE": "A",\
              "URA": "U", "DA": "A", "DT": "T", "DC": "C", "DG": "G"}

AMINOACIDS = {"A": "ALA", "R": "ARG", "D": "ASP", "C": "CYS",\
              "E": "GLU", "G": "GLY", "Q": "GLN", "H": "HIS", \
              "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET",\
              "F": "PHE", "P": "PRO", "S": "SER", "T": "THR",\
              "W": "TRP", "Y": "TYR", "V": "VAL", "N": "ASN" }

NUCLEOTIDES = {"C": "CYT", "T": "THY",\
               "G": "GUA", "A": "ADE",\
               "U": "URA"} #check if we need to add modified nts names

AANAMES =    ["ALA", "ARG", "ASP", "ASN", "CYS",\
              "GLU", "GLY", "GLN", "HIS", \
              "ILE", "LEU", "LYS", "MET",\
              "PHE", "PRO", "SER", "THR",\
              "TRP", "TYR", "VAL"]

NTNAMES = ["CYT","THY", "GUA","ADE", "URA"]

#"DA" and "CA" are both nts and ligands; removed from ligands

LIGANDS = ['M3L', 'SO4', 'HEM', 'DAL', 'MLE', 'MVA', 'BMT', 'ABA', 'SAR', 'NH2', 'NA', '00K', 'TYS', 'IGU', 'MSE', 'FME', 'ZN', 'SCH', 'CA', 'THP', 'EFC', 'GGL', 'CYW', 'SEP', 'TPO', 'OCT', 'DAR', 'DGL', 'DIL', 'DXX', 'MGG', 'CSS', 'T42', 'U1', 'NLE', 'DBB', 'MG', 'GDP', 'ACE', 'HG', 'TRI', 'GTP', 'GPL', 'MN', 'PCA', 'NAG', 'BMA', 'CD', '1PG', 'BZT', 'DPP', 'LLP', 'GLS', 'NO3', 'ETA', 'MPD', 'CL', 'MRD', 'PI8', 'BRU', 'FVA', 'DLE', 'DVA', 'ACY', 'DPN', 'HMR', 'NAD', 'ADC', 'IH1', '1AP', 'HYP', 'G49', 'DM1', 'ATP', 'CGU', 'CXM', 'MAN', 'XYS', 'IVA', 'DFO', 'NME', 'DMF', 'NVA', 'DCY', 'DHI', 'DTR', 'BR', 'CYN', 'ACT', 'CSO', 'FUC', 'T3P', 'U34', 'TPQ', 'CU', 'CRO', 'G4S', 'DGS', 'DMT', 'X', 'LGP', 'CSD', 'COA', '5CM', 'CUL', 'FGL', 'NDG', 'PTR', 'ANP', 'DOC', 'DCT', 'MLU', 'OMZ', 'GHP', 'OMY', '3FG', 'BGC', 'RER', 'LAC', 'ALZ', 'PHD', 'C45', 'PO4', 'DPM', 'PTH', 'DIP', 'DSN', 'AP2', '1ZV', 'PI9', 'CSW', '612', 'NI', 'G3P', 'NO', 'SNC', 'MOH', 'AAB', 'BCT', 'XYP', 'BPP', 'AS', 'DM3', 'YRR', 'TMD', 'MSA', 'LML', 'DGN', 'DPR', 'DAS', 'MHL', 'T19', '0IV', '00L', 'PM3', 'U31', 'GAL', 'AAL', 'IH2', 'ARB', 'PLP', 'NGU', '64T', 'SC', 'SOR', 'GLV', 'PRL', 'DTY', 'CY3', 'FMT', 'DUR', '5AT', 'TRP', 'UMP', 'MTX', 'CO', 'SAH', '7MG', 'PSU', '5MU', '1MA', '5MC', 'TYY', 'HY1', 'PEO', 'PEA', 'GOL', 'TIH', 'IGL', 'OIC', 'CF0', 'EDO', 'HSY', 'DFI', 'MDL', 'PLA', 'STA', 'EHN', 'ACB', 'MFD', 'FGA', 'MDH', '0KV', 'HTR', 'CB3', 'FE', 'PME', 'CLD', 'EDC', 'CLA', 'PQN', 'SF4', 'QWE', 'MYR', 'IOD', 'OES', 'AHO', 'CSP', 'MNV', 'TBM', 'NAP', 'SBI', 'CS', 'BT2', 'PYR', 'TRA', 'GSE', 'BLA', 'MEN', 'CYC', 'BO4', 'GLC', 'SVA', 'ARG', 'CSI', 'PHA', 'APA', 'CYQ', 'CSE', 'FE2', 'FCO', 'H2S', 'CPT', '3DR', 'FTY', 'SEB', '6HG', '6HT', '6HA', '6HC', 'RTY', 'KCX', '3PG', 'A23', 'K', 'ADP', 'ORN', 'NET', 'CYG', 'UFP', 'TMF', 'PI7', 'NPH', 'BUC', 'AYA', '5IU', 'IAS', 'ALC', 'MPQ', 'TBG', 'BAI', 'CY1', 'DHD', 'ANS', 'AR7', '0QE', 'MGD', '4MO', 'GCA', 'BAM', 'G36', 'NT', 'NGM', 'BM5', 'IIC', 'HIC', 'FOE', 'CTH', 'POM', 'MAA', 'PXZ', 'H5M', 'BNZ', 'P', 'IMC', 'C38', 'DNP', 'PI1', 'TOL', 'TYQ', 'R', 'CIT', 'CME', 'U10', 'PEE', 'BOG', 'FES', 'FAR', 'CRF', 'PYX', '6MA', 'GSS', 'DIY', 'RE', '34H', 'PRJ', 'OAR', 'ATD', 'HAQ', 'SEL', 'RET', 'CMC', 'LMR', 'DAB', 'C49', 'T32', 'ARC', 'SGA', 'NCO', 'SIN', 'ALN', 'C31', 'MMT', 'AIB', 'PHL', 'CCN', 'G47', 'IML', 'BA1', 'SLZ', 'ODA', 'EPE', 'IPA', 'HED', 'DDU', 'LAA', 'MCY', '0A9', 'GLA', 'PT', 'HTO', 'BB9', 'BB6', 'BB7', 'BB8', 'MH6', 'TSP', 'TCP', '00N', 'CBR', 'C6C', 'IBZ', 'GPE', 'GC4', 'ASG', 'GCU', 'A40', 'MAR', '4SU', 'H2U', 'QUO', 'G7M', 'AMP', 'AMO', 'CDE', 'PI5', 'Y', 'TRS', 'MNL', 'GSR', 'ESI', 'SPM', 'IYR', 'LTA', 'TRN', 'AGL', 'HMC', 'BAB', 'PI4', 'PHI', 'SMC', 'NOR', 'CAP', '00P', 'SET', 'CAS', 'DIX', 'QIL', 'EOH', 'GR4', 'BT1', 'DLY', 'S02', 'SR', 'THX', '6OG', 'PEB', 'PUB', 'AIR', 'BDP', 'CEG', 'GLZ', 'XPR', 'AG2', 'SRA', 'BM1', 'AZL', '0G6', 'BHD', 'IH3', 'BIN', 'CFF', 'DSE', 'GMA', 'RIP', 'IUM', 'HPE', 'POL', 'OCS', 'DGP', 'SGN', 'IDS', 'GAM', 'EEE', 'BT3', 'SME', 'BAH', 'G31', '0Z6', 'CAC', 'SBL', '1PA', 'CSB', 'SAC', 'E', 'MIS', 'ACN', 'PI6', 'C5C', 'DFT', '1MG', 'FLT', 'AEA', 'VAD', 'MTY', 'OH', 'PRN', 'OCY', '0ZE', 'MIA', 'CYS', 'GNP', 'A38', 'DPG', '00R', 'SNN', 'PSE', 'L2O', '0IT', 'TYI', 'PED', 'ILX', 'TRX', 'CSX', 'CPL', 'U5P', 'C5P', '5GP', 'LA', 'ALY', 'MES', '2DF', 'AC1', '5PC', 'PDU', 'HIS', 'BAL', 'CMO', 'FTR', '13P', 'CRG', 'AD2', 'PPU', 'CSR', 'U', 'PDI', 'PG4', 'PAR', 'MLY', 'V4O', 'T76', 'DBH', 'ADE', 'CDD', 'IAP', '0AF', 'HDZ', '669', 'B2A', '5BU', 'ACO', 'POP', 'CYD', 'SO2', 'P5P', 'VPR', 'POA', 'PLM', 'CXS', 'S6G', 'SMA', 'CH', 'LYX', 'CDI', 'ENX', '2MG', 'M2G', 'OMC', 'OMG', 'YG', 'MAD', 'UDP', 'OMT', 'O', '1PE', 'PEG', 'PLG', 'DIV', 'DDX', 'UCN', 'PIB', 'HPD', '0HQ', 'Y1', 'SAM', '1OL', 'TPP', '15P', 'BGM', '8OG', 'HOQ', 'PHE', 'PEP', 'DCS', 'BME', 'KI2', 'CO3', 'THG', 'DTU', 'DTT', '4OX', 'CIR', 'PFF', 'MEL', 'DRZ', 'ANG', 'MTR', 'CMT', 'HDS', 'PPI', 'DHA', 'TEE', 'HEU', 'FHU', 'QPH', 'ABE', 'NTC', 'DOF', 'GD3', 'ATL', 'ST8', 'SGM', '2AS', 'MAL', 'ADW', 'HDU', 'IGN', 'TYB', 'G6S', 'NGS', 'SPS', 'MAG', 'NIY', 'FAD', 'C24', 'ABF', 'XYL', 'MCL', 'LEF', 'HSM', '5CG', 'GPN', 'TPN', 'APN', 'A66', 'T66', 'C66', 'CPN', 'B2V', 'GBX', 'NHE', 'MHS', 'AGM', 'MGN', 'GL3', 'F43', 'TP7', 'COM', '6MO', 'CDL', 'HQO', 'MPT', 'CLG', 'CLH', 'DMS', 'FDC', 'LDA', 'CDB', 'DSB', 'PDD', 'THC', 'CYF', '2PE', 'IIL', '790', 'ESD', 'YT3', 'DOE', 'OLA', 'BEO', 'UNX', 'SCY', '1BO', 'SUJ', 'CNT', 'GDX', 'DLP', 'DNR', 'MC1', 'SUI', 'SER', 'THR', 'MRG', 'ATM', 'MN8', 'MN1', 'CHD', 'CMR', 'RMP', 'SMP', '3ME', 'HDT', '4IN', '5ZA', '16A', 'A2M', 'FMN', 'THJ', 'ASP', 'CB2', 'HTI', 'GCK', 'MOE', 'SPK', 'VOL', 'TRW', 'URA', 'LAE', 'SO3', '01B', 'AVC', 'NVP', 'CR9', 'DMK', 'DMH', 'CFT', '3OH', 'ALF', 'AF3', '233', 'SQ', 'BRG', 'IMD', 'MLT', 'PST', 'TLP', 'TLB', 'NBS', 'HIP', 'NYC', 'MDP', 'UMS', 'NIT', 'UBB', 'GNG', 'CAR', 'GCP', 'TYN', 'BP4', 'OXI', 'RG1', 'BCL', 'BEN', 'CG2', 'STU', 'DP5', 'TRO', 'HEC', 'TFB', 'PQQ', '2DA', 'TRQ', 'XE', 'NC1', 'AHB', 'SPT', 'YCM', 'U3P', 'ADN', 'U05', 'KIV', 'CY4', '162', 'SMF', 'CDM', '5CS', '1PR', 'PTY', 'MC', 'AKG', 'RHO', 'C02', 'MTL', 'CBI', '696', 'T2S', 'C2S', 'G2S', 'G4P', 'TFT', 'CR2', 'GA4', 'ATG', 'N76', 'Z', 'CRQ', 'ACP', '163', 'CCS', 'ETC', '4ND', 'LPS', 'ABS', 'DCP', 'TTP', 'HDY', 'OMU', 'CCC', 'TGP', 'TTG', 'TTC', 'MLZ', '168', 'HSO', 'ESC', 'CSZ', 'RAM', 'MXZ', 'DGC', 'NGK', 'HAL', 'A3P', 'MHO', 'CDA', 'BOR', 'FGP', 'SUC', 'TYT', 'I59', 'A2G', 'KIR', 'IMP', 'L86', 'ADZ', '1ZN', 'DAM', 'FCL', 'MOA', 'SCN', 'PMP', 'CHI', 'CSY', 'QUA', 'DBU', 'TS9', 'FSC', 'SPP', 'PLD', 'UPL', 'N20', 'BOC', 'LOV', 'TLA', 'BNG', 'MMC', 'SMT', 'CNC', 'AA4', 'C8E', 'RVP', 'ABR', 'PG1', 'C36', 'IU', 'NDP', 'MN7', 'MN2', 'PP3', '143', 'N41', 'FSN', 'M1G', 'SSU', 'ASU', 'BLE', 'CAA', 'ORP', 'LAR', 'S1H', 'FUL', 'HAI', 'YYG', 'NEP', 'YCP', 'TLX', 'TLC', 'R11', 'PGO', 'HSA', 'TTN', 'U8U', 'T6A', 'INR', 'BCS', '2MR', 'THF', 'FTC', '2AR', 'UPF', 'MAK', 'FSO', 'F3S', '3PA', '4BA', 'TYR', 'PHO', 'NO2', 'SNP', 'EFZ', 'CSA', 'ATO', 'PHQ', 'CP2', 'DII', 'UFR', '8AD', 'MBZ', 'CG1', 'BJI', 'KDO', 'GMH', 'GPH', 'GP4', 'GP1', 'LIL', 'AAE', 'LIM', 'EA2', 'FCI', 'PGP', 'LU', 'OXL', 'TAR', '132', '112', 'BE2', 'DBY', 'DCE', '1Z0', 'MPR', '130', 'TZE', 'EDA', 'BBS', 'CPI', 'DOA', 'TPR', 'TP2', 'LYN', 'ORD', '2TL', 'D4P', 'ALO', 'CHP', 'FA7', 'RAL', '3AT', '3PO', 'FPT', 'MSO', 'PNI', '4IP', 'NRI', 'A47', 'CSH', 'G1P', '0DC', '0DG', 'NEA', 'ARO', 'CDP', 'P2U', 'B2F', 'T87', 'HEO', 'SIA', 'NGA', 'TEN', 'PEL', '121', 'ASA', 'DEL', 'PN2', 'HPH', 'OAS', 'SCL', 'WO2', 'NW1', 'DCL', 'PU', 'BCB', 'BPB', 'MQ9', 'MST', 'NS5', 'T23', '0DT', 'TP3', 'IG', 'IC', 'S11', 'CMG', 'PYC', '9AC', 'UAG', 'API', 'CCY', 'DMQ', 'MTG', 'LNK', 'I84', 'BMZ', '5HT', 'HAE', 'ROS', '0ZI', '2MA', '120', 'PBF', 'ADQ', 'GER', 'PYP', 'AKB', 'ALS', 'GFL', 'TAF', 'CFL', 'FRU', 'GTT', 'TP4', '1LU', 'PYY', 'EPC', 'PSA', 'LYK', 'ATH', 'MME', '3AH', '2AN', '0QS', 'NP3', 'DLS', 'OSE', 'DDZ', 'CGA', '2GP', 'LPA', 'TPL', 'ZAL', 'MMO', 'PCD', 'ORO', 'G', 'CSU', '2PP', '0MG', '70U', '12A', '2MU', '3Q5', 'UNL', '64U', 'G6P', '5ID', 'UMQ', '8PE', 'CN6', '9PE', '7PH', '6PH', 'FSA', 'RB', 'QQ2', 'G2M', 'COD', 'DM0', 'MCM', 'PUK', 'NFA', 'IEY', 'MZZ', 'BTB', 'AGT', '2EG', 'LYZ', 'DBV', 'O2C', 'C1M', 'FDA', 'PE3', 'LPD', 'TLN', 'LCG', 'J2I', 'HBQ', 'JG3', 'VME', 'LYS', 'NRQ', 'CH6', 'EYG', 'UIR', 'HLU', 'MNB', 'IZA', 'RIA', 'DTP', 'D94', 'DDE', 'APR', 'SO1', 'MFT', 'SCC', 'IPE', 'OBN', 'HY3', 'GLU', 'XMP', 'EQO', 'MPO', '6CW', '011', 'B7C', 'KYN', 'BU1', 'PXY', 'AGQ', 'MA7', 'XL3', 'PPX', '1EM', 'ARX', 'HRG', 'PGE', 'MA5', 'ASB', 'T04', 'AZI', 'CR0', 'NH4', 'P1L', 'EBP', 'EFS', '7DE', 'CSF', 'PRF', 'O12', 'LMT', 'PB', 'PRP', 'OBS', 'N6R', 'N6S', 'BAP', 'RUS', 'UIB', 'IT1', 'PZA', 'AZA', '2BT', 'ABU', 'ACA', 'NAL', 'CP', 'PL9', 'OEC', 'BCR', 'DGD', 'LHG', 'SQD', 'LMG', 'FHO', 'PVE', 'BI2', 'BU3', 'COK', 'A5N', 'DG2', 'C42', 'G38', 'A43', 'NYM', '4SP', 'TG1', 'DYG', 'G3B', 'CBV', 'FMR', 'GLY', 'ALA', '10M', 'DRR', 'DPJ', 'S4C', 'UQ1', 'E1X', 'M0E', 'OHI', 'OXY', '3DM', 'D2M', 'KIA', 'CMK', '5FC', 'AHP', 'DIO', 'DOR', 'MLI', 'NCD', '6SA', 'AYG', 'CQP', '2AD', 'TRE', 'URN', 'HPP', 'TBS', 'SX7', 'BI3', 'PE4', 'OGA', '2NC', 'MY5', 'PBT', 'THM', 'R1P', '376', 'AE3', 'CO2', 'D33', 'AG', 'AX4', 'D1D', '23F', 'FLC', 'DXD', 'SHA', 'TQQ', '2DT', 'EOD', 'ARS', 'RAJ', '0GR', 'ANL', 'F42', 'LCC', '897', '2PO', 'HDD', 'F50', 'D3T', 'DAP', 'UN9', '2PR', 'LET', 'CMY', 'BA', 'DGT', 'BZG', 'NVB', 'YMP', '6IA', 'TAN', '49U', 'GVJ', 'IPH', 'MEA', 'ES4', 'V1A', 'HEZ', 'PGR', 'AZZ', 'OAK', 'SMM', 'KGC', 'D91', 'CYI', 'KSS', '11U', 'BBL', 'QMM', '16G', 'ALT', 'C2F', 'XDR', 'DDG', '3MY', 'T55', 'HEA', 'TGL', 'PGV', 'CUA', 'PEK', 'PSC', 'DMU', '4AF', 'CTP', 'GIG', 'PR9', 'DM8', '4CY', 'XCY', 'TAM', 'SN3', 'NDB', 'C3M', 'PR4', 'IAM', 'FOX', 'DNG', 'DNE', 'DNM', 'NLO', 'XP9', 'MAE', 'TES', 'B0D', 'B9D', 'ASJ', 'DDS', 'MDO', 'GFF', 'MLL', 'ORQ', 'XYG', 'S9L', 'MTU', 'B3T', 'XCP', 'XPC', 'B3E', 'N2C', 'NCY', 'QUI', 'DZ4', 'TTD', 'FO1', '5PY', 'HC4', 'DEU', 'PMI', 'DTV', 'FXP', 'GR1', 'DTD', 'GAH', 'MH4', 'SPX', 'O2G', 'AF', 'DAH', 'TY2', 'UR3', 'MYL', '5AA', '2OP', 'PO2', 'CM0', '3Q2', 'GVL', '4PP', 'DGI', '32T', 'M1P', '24S', 'URE', '173', 'ARF', 'ROC', 'A3D', '0C', '0U', '0G', 'ADI', 'MUB', 'FPP', 'OTT', 'MP8', 'UPS', 'IHP', 'TYZ', 'LEI', 'ZZJ', '8AG', 'OAD', 'LI', '3QN', '3Q0', 'F6P', 'D3P', 'R1A', 'YB2', 'CTG', 'SAP', 'CXP', 'WFP', 'YNM', 'JOZ', '0AH', 'SAL', 'G3H', 'CGN', 'FMX', 'UQ', 'SRT', 'LED', 'BEZ', 'FMG', 'OZ2', 'NEH', 'BEF', 'XFE', 'PP1', 'AX5', 'C2B', 'PG6', 'A5M', '8LI', '66G', 'CE6', 'SMY', 'KBE', 'UAL', 'MYN', 'NU5', 'LEA', 'HT', 'UC1', 'MET', 'OSM', 'ILG', 'FOT', '211', 'OTG', 'LK2', 'TB6', '5PI', 'EPS', 'G3D', 'M2T', 'SEC', 'M04', 'P6G', 'A1R', '459', 'MP5', 'PA9', 'LNT', 'SHT', 'FRD', 'HEF', 'UDC', 'US1', 'F', 'GPD', '2KT', 'TMP', 'DHF', 'MUA', '23T', 'WO4', 'MMA', 'IQS', 'GYC', 'SPD', 'VO4', '382', 'DOL', 'MHW', 'MHU', 'MHV', '004', 'MHT', 'BFD', 'P4C', 'G12', 'NPO', 'V25', 'VIR', 'C12', 'NMN', '18S', 'LLO', 'L9L', 'PZP', 'CCR', 'C2A', 'PG0', 'XP8', '17M', 'LE1', 'C2E', 'DRQ', 'DRW', 'MD1', 'AGA', '3PH', 'DRK', 'ECX', 'CEJ', 'A5P', 'R4A', 'LEU', '02G', '03Y', '1ZZ', 'TIG', '3AZ', 'BG6', 'LA2', '2HP', '894', 'B3M', 'BIL', 'B3L', 'XGR', 'XCR', 'XTR', 'XAR', 'TSS', '2MD', 'MO', 'UMA', '18O', 'PYO', 'DCC', 'DAO', 'PIL', '3GC', 'N8E', 'FH7', 'PF5', '12Q', 'B3K', 'B3D', '8AN', 'OXM', 'TAD', '23U', 'AIN', 'R69', 'XX5', 'AX2', '4FB', 'SM', 'JK2', 'UNC', 'NNR', 'DAD', 'SE', '5HU', 'CMH', 'DXC', 'PSH', 'JZH', 'ROL', 'XXY', 'F25', '3CO', 'ZHH', 'ZHZ', 'BDR', '2FG', 'T5S', 'GWI', '6DE', 'ANM', 'GMP', 'MG1', 'MBN', 'SLU', 'NYS', 'MQ7', 'NS1', 'DEP', 'NFO', 'FRF', 'SFG', 'OXG', 'TA2', 'XYQ', 'NCA', 'MNU', 'MK8', '4LM', 'CRK', 'G95', 'MXE', 'Z3E', 'ECQ', 'KCQ', '0CS', 'DG3', 'GAI', 'KOR', 'K1R', 'H1L', 'Z57', 'DHK', 'RC8', 'TTI', 'UPG', 'TOX', 'HOX', 'PPV', 'TM1', 'MAU', 'TSR', 'GYS', 'CML', 'CAF', '976', 'P03', 'PER', 'CR7', 'ABG', 'G3A', 'OTY', 'CPS', 'SAT', 'N2G', 'CIC', 'A1P', 'AHR', '4MF', 'S06', '4HT', 'ONE', 'LIV', 'F2A', 'IRI', 'DHB', 'BTN', 'A0A', 'EXC', 'MY2', 'JAS', '6PL', '6UL', '6JZ', 'CRU', 'BIF', 'VLM', '4SC', 'DOD', 'B83', 'TSE', 'C', 'ABY', 'ZYX', 'AME', 'LYR', 'ICT', 'CLB', 'SCA', 'P63', 'M03', '896', '53U', 'GVP', 'BRN', '051', 'DBZ', 'KST', 'PSZ', 'FSE', 'TDP', 'P33', '2SM', 'FAB', 'TDM', 'CVM', 'XIX', 'OAA', 'SCX', 'N5C', 'N5M', 'L2G', '2ST', 'CYR', 'PAQ', 'SM9', 'GVH', 'B99', 'X2L', '887', 'M07', 'SAI', '8MG', 'PI0', 'PVX', 'GRB', 'CBS', 'XM1', '0AR', 'UN1', 'ZGU', 'ZBC', 'ZCY', '2FI', 'MOO', 'DPQ', '0MO', 'IPN', 'LKM', 'HHK', 'GSH', 'UQ2', 'MGF', 'A5L', '2GT', 'G98', 'SUN', 'PA1', 'GCN', 'FTT', 'DPO', 'EAP', 'DDQ', '00Q', 'GNE', 'B8D', 'GLF', '10C', 'LCA', 'ESH', '37U', 'CAM', 'AH0', 'FOA', 'MQA', 'EAB', 'KPI', 'ZZI', 'CLE', 'PRO', 'D5M', 'SPE', 'B33', 'OSC', 'ZDU', '7GU', 'TCA', 'PG5', 'OXC', 'PX1', 'OLN', 'GSP', 'GF4', '179', 'HQA', 'PR', 'VR0', 'M2L', 'PVW', 'N2O', 'XXX', '2TS', 'HSE', 'IBN', '4PC', '23D', 'SFU', 'CR8', 'IEL', '4OC', 'MA6', 'XX1', 'DYS', 'CFY', '57D', 'EUG', 'VLL', 'VDL', 'LNG', '3SP', 'ETX', 'G25', '0ZS', 'SSJ', '3CS', 'P19', 'LVN', 'SWI', 'CP1', 'TP1', 'S13', 'GWB', 'GGH', 'T3', 'MTA', 'RA8', 'TNV', 'TA4', 'HL2', 'OTH', 'HF2', 'GUN', 'LLY', 'CHG', 'BUA', 'PO3', 'KAG', 'JZV', 'SGB', 'HSX', 'CZP', '3TY', 'BZZ', 'VIV', 'TB0', 'IIP', 'DMP', 'YTP', 'T12', 'G96', 'APK', 'US2', '925', 'SCJ', 'LM2', 'CRX', 'MMZ', 'TYD', 'C0T', 'NRP', 'ACM', 'BNO', 'URC', 'H1S', 'PA5', 'TXZ', '06C', 'GLP', 'PRV', 'PRQ', 'CTZ', 'G16', '2CE', 'ME6', '40G', '40C', '40A', '40T', 'FP9', 'FTE', 'CPB', 'AIS', '6MN', '69P', 'OCE', 'GHG', 'MX3', 'PXP', 'OI7', 'B3S', 'B3Y', 'B3X', 'B3A', 'VX6', 'PHZ', 'RIO', 'DT6', 'OPE', 'UNK', 'MLC', 'P4G', 'CZB', 'DA2', 'OLU', 'QUE', 'HFS', 'ACH', '608', '22S', 'SN6', 'PZE', '2HG', 'S91', 'SDG', 'ANY', 'PLC', 'TY5', 'RE0', 'ABN', 'EXB', 'CFZ', 'AF2', 'UFT', 'CQU', 'BK2', 'C46', 'AAR', 'NLK', 'SCS', '3CM', '6SC', 'TAY', 'CW1', 'CLJ', 'T5O', '3SA', 'OPR', 'MTN', 'BLM', 'FIC', '4ST', 'AUK', 'PS9', 'DGB', 'RC7', 'CA1', '1P1', 'SC2', 'HEX', 'UQ8', 'KAB', 'B3Q', 'TBU', 'LL2', '5B2', 'TTQ', 'B22', 'ACR', 'ITT', 'BET', '689', 'LK3', 'OHX', 'HL4', '247', 'CZA', 'X4A', 'F55', 'SXE', 'R5A', 'R5B', 'AGS', '295', '51U', 'XSN', 'GLH', 'BAT', 'PAU', 'PFU', '2RL', 'LFR', 'PLR', 'GDD', 'DW2', 'TB', 'VDM', 'OPI', 'RUG', 'FLD', 'MLA', '3FP', 'FBP', 'R1F', '55E', 'NJQ', 'W', 'U33', 'BCN', 'P47', 'HBB', 'BCD', 'Z77', 'MI1', 'MXC', 'LYP', 'RUB', 'HS2', '195', 'PE5', 'P42', 'A03', 'CL1', 'RFZ', 'L03', 'PXF', 'BF4', 'BF2', 'TYE', 'A93', 'UTP', 'TT', 'XLS', 'DW1', 'L3O', '7PE', 'XJS', 'OCQ', 'J3Z', '1DP', 'SCE', '193', 'GVD', 'L3G', 'PQ9', 'MGE', 'AZH', 'STH', 'GWE', '4DP', 'FZN', '5FU', 'PCI', 'NCS', 'CS1', 'FAG', 'GVI', 'PGA', '32U', 'OZT', 'WO3', '63G', 'TL', 'BIB', 'TFE', 'DKA', 'HOL', 'CH7', 'I3A', '0AK', 'I5S', 'TAU', 'APC', 'EU', 'PDC', 'IBF', 'R7A', 'F68', 'DCZ', 'HFA', 'LIF', 'DKY', 'WQP', 'TDG', 'JDT', 'IFM', 'DKX', 'J30', 'BI4', 'CEB', 'PEF', 'LPP', 'MUL', 'HOA', 'H52', 'HBP', 'IAC', 'M4C', '3MU', '3AU', 'CN5', 'CN3', 'NBT', '2ML', 'KAP', 'DCQ', '6MZ', 'SAY', 'SN0', 'HNF', 'HS8', 'ML3', 'RE3', 'TYL', 'MSU', 'ALV', 'PIH', 'SYS', 'DSG', 'DTH', 'GYV', 'PC3', 'BPE', 'CYT', 'BI1', 'QPS', 'DTQ', 'IPG', '1SC', 'C43', 'G48', 'A44', 'U36', 'C4R', 'RH', 'KPN', '3CF', 'TC1', 'K66', 'AM9', '3CN', 'BHO', 'PUT', 'M8M', 'AU', 'MG8', '0MA', '157', 'G11', '0GJ', 'B30', '00C', '4E6', 'HTG', 'ID3', 'USM', 'QRW', 'NIA', '3DA', 'MIN', 'JZA', 'TRT', '5PG', 'M12', 'FOR', 'T2T', 'CMP', '1TB', 'P22', 'NIM', '22U', 'GVE', 'TFO', '1AC', 'HFP', 'T4S', 'ZHP', 'ZAD', 'ZTH', 'HCI', 'TYM', 'FIL', 'A96', '3TL', '2AU', 'AQC', 'PDS', 'VAL', 'LAT', '3RC', '0BN', 'FUM', '3TR', 'K12', 'DFP', 'J60', 'AHD', 'IPR', 'ZS0', 'FU4', 'DMG', '4HL', '19U', 'GMB', 'SOS', 'M1T', '985', 'SMR', 'UIP', 'SUR', 'MED', 'BO3', 'PH2', '3Q6', 'UAR', 'C5M', 'BPI', 'DX7', 'BDH', 'GAB', 'OXX', 'TBZ', 'FEO', 'MP6', 'LME', 'PZR', 'LGY', 'CQW', 'TAC', 'IMN', 'OHA', 'P44', 'HNC', '501', 'ARE', 'PP2', 'AA5', '13U', 'M05', 'PPN', 'XX7', 'JIM', 'RHM', 'GLN', 'KBM', 'A3A', 'KZM', 'STR', 'SCR', '380', 'DCR', 'JK3', '50U', 'IKR', '538', 'GK5', 'PE8', 'PVA', 'AP7', '01W', '2NT', '4PH', 'PAE', 'JW5', 'MFX', '4OH', 'M38', 'EPU', 'HBI', 'CQA', 'AFB', '1PS', 'CY7', 'C75', 'NNN', 'GLK', 'RA4', 'VDN', 'PDO', 'ZUK', '5SC', 'MHR', 'DIF', 'XX6', 'UC3', 'TOF', 'HCX', '45U', 'ITI', 'M01', 'CFD', 'G35', 'I08', 'KSK', 'UC2', 'CK4', '4HY', 'A02', 'VLT', '0FL', 'AHH', 'DA6', 'OTD', 'MIR', 'H18', 'APJ', 'HMT', '771', 'GRG', '1MM', 'Z82', 'DSD', '63H', '2MT', 'F2F', 'M5M', '3CB', 'HI6', 'DCF', 'PFX', 'TA6', 'YOF', 'ETM', 'AT3', 'GFA', 'FFD', 'R7U', 'FNM', 'FPA', '91U', 'DMA', 'AKZ', 'HZP', 'ONL', 'P3Q', 'SVZ', 'YB', '5AX', 'BOE', 'BTR', 'GTH', 'HBH', '3Q4', 'AT2', '24U', 'RNG', 'NDU', '6PB', 'EA1', 'TSC', '19B', 'T39', '3BP', 'MXF', 'EIT', 'M25', 'PIA', 'PLZ', '5OH', 'SIC', 'N12', 'GCS', 'CJO', 'PEC', 'CZ2', 'TAS', 'URI', 'DRS', 'RP5', 'MCS', 'EPL', '857', 'PD8', '3TD', 'JK1', 'DM6', 'I26', 'FFQ', 'UD1', 'I25', '5B3', 'L17', 'PR7', 'BH2', 'LSP', 'ESB', 'K68', '296', 'B72', 'C1X', 'DPF', 'TBF', 'C5X', 'XIZ', 'GIC', '166', 'TED', 'DMR', 'XGL', 'XTL', 'XAL', 'XCL', 'FS2', '25E', 'POG', 'AMZ', 'IDC', 'M29', 'AD6', 'F89', 'CHX', 'MDJ', 'TFA', '4CF', 'SIG', 'AMY', 'XUG', 'NIO', 'T8N', 'R94', '126', '255', 'LAV', 'G1C', 'RHD', 'LI8', 'DPL', '64P', '8FG', 'H12', 'SEN', 'WR2', 'OEG', 'LCK', 'SLB', 'TYX', '0AD', 'RSP', 'KDA', 'NZH', '5UA', 'RXV', 'ZBZ', 'MKC', 'TSH', '12U', 'GK6', 'FOL', 'PQA', 'FDP', 'IM9', 'VIB', 'NAI', 'VGD', '743', '2OT', '16O', 'DL7', 'H8H', 'KZL', 'MGL', 'SGC', 'SCF', 'IHS', '3Q3', 'HPX', '24O', 'LYT', 'PNZ', 'KYQ', 'KSL', 'SOC', '176', 'OSB', '33U', 'CHS', 'HPA', 'CAQ', 'MDV', 'M6P', 'FER', 'UVX', 'SB2', 'TXQ', 'CE1', 'NFC', '2HE', 'CU1', 'M18', '2HF', 'EHB', 'BBC', 'GG7', 'L0I', 'KOT', 'NBG', '055', '5CF', '03R', 'LAL', 'HCS', 'MRY', 'GCO', 'COI', 'PZM', 'SSS', 'N6G', 'KZI', 'NF2', 'GME', 'FCA', 'MOI', 'FRV', 'A5A', 'HS9', '1P5', 'TDY', 'CDW', 'G27', '3SL', 'HBA', 'LDP', 'P5B', 'HY0', 'XIY', 'HOB', 'PSM', 'Z15', 'Z16', 'B3P', '7PA', 'ELA', '44U', 'T3O', '12E', '6MI', 'DK4', 'DHL', 'MFN', 'H4Z', 'MPG', 'O8M', 'PTN', 'TNK', 'ER3', 'NMS', 'THO', 'I19', 'M3R', 'LK4', 'GOA', 'GVG', 'R1B', 'LDH', 'HSL', 'P16', '7DA', 'R68', 'XXA', 'HXZ', 'ZBA', 'XG4', 'AFC', '22O', 'NMC', 'PR3', 'X39', 'CSL', 'DAI', 'LAG', 'GOX', 'CK7', 'AHY', 'YJD', 'IN5', 'P2T', 'AZO', 'SAV', 'F6R', 'ULA', 'TTS', 'GVN', 'C85', 'CTC', 'C8M', '475', 'LLZ', 'Y27', 'UDA', 'CET', '1AB', '2FE', 'ZRL', 'CLV', 'B15', 'SVR', 'VIA', 'MEQ', 'PPW', '165', 'LLL', 'PEU', 'AUC', 'BA2', 'ZBU', 'AEI', 'HMH', 'DKI', 'KSR', 'TCQ', 'BI8', 'OHN', 'KIM', 'C6P', 'MX5', 'L1G', 'PMH', '0ZT', 'CKG', 'TPS', 'TOH', '2C2', 'NA8', 'TH6', 'PRS', 'X2K', '12B', '2BA', 'T49', 'GMU', 'OMH', 'HMG', 'NH3', 'E64', 'CG3', 'LWY', 'M8E', 'GPR', 'CRS', 'SVV', 'SVW', '20X', 'RXA', 'TKL', 'NIC', 'PPY', 'SKE', 'SMH', 'SCQ', '21U', 'VX2', 'BIO', 'C94', 'GBR', 'ZAB', 'HXB', 'PXL', 'ICX', 'HRM', '2LG', '3GL', 'NO1', 'H16', '348', 'R9A', 'D11', 'US4', 'JJL', 'AX6', '13T', 'SVY', 'X1N', 'IZO', 'H20', 'OHS', 'C4X', '0AS', '46U', 'AX1', 'MHF', 'CSC', 'GF2', '194', 'TTM', 'L20', 'BFC', 'L0G', '9BD', '914', 'MXS', 'NFR', 'RXC', 'GFC', '900', 'AAC', '29U', 'AHZ', 'XS2', '0ZC', 'AMU', '464', 'EGC', 'FUA', 'D15', 'N10', 'OMN', 'TON', 'SDP', 'SCZ', 'CNA', 'MCO', 'C4S', 'S4G', 'S4A', 'PC0', 'URX', 'EEM', 'SWE', 'EDN', '17T', 'DL8', '627', 'PQ1', '9DG', 'FA2', '177', 'SSA', '0G7', 'T41', 'ILE', 'MDK', '575', 'HIA', 'FSX', 'SBY', 'PSW', '2PM', '16F', '4DB', 'DX2', 'LNC', 'MDU', 'EDD', 'PFI', '723', 'XUA', 'STE', 'SBG', '3XH', 'EST', 'ALQ', 'GLJ', 'TYJ', '3Q1', '5IQ', 'NMM', 'LMU', '1P2', 'SVX', 'GS', 'P55', 'IQB', 'TB9', 'T15', 'NYB', 'CPC', 'HQU', 'SX', 'FOU', '4AX', 'NEN', 'P04', 'IN1', 'KDP', 'BZD', 'R5P', 'ANZ', 'TIS', 'T0I', 'NOJ', 'NMY', 'RIB', 'M08', 'DPV', 'TR9', '167', 'TCY', 'DTB', 'X1P', '1AE', 'AA1', 'NBQ', 'HT5', 'BYH', 'AN0', '3PB', '125', '01I', 'DX6', '9AT', 'ASN', '5CB', 'SA3', 'AMG', '26B', 'CZO', 'FX1', 'B1L', 'IP8', 'ZP4', 'NA9', 'AOR', 'SGR', 'SGX', 'DZM', '1TQ', 'DVT', 'HLO', 'MMR', 'A8M', 'GR0', 'LEW', 'CIE', 'H1N', 'ET0', 'BGG', '12M', 'MFR', 'PRK', 'QBT', 'GYT', 'US5', '4BF', 'HAO', '4DG', 'PPF', '1MS', 'AP8', 'X4Z', '2BU', 'PCF', 'Q6W', '5BD', 'HT1', 'MY1', 'OGX', 'GA9', '2BR', 'XYA', '0PY', 'PIE', 'DRA', 'TY8', 'TY9', 'ME0', '567', 'UNB', '16P', 'AR4', 'HIQ', '706', 'JN5', 'EFP', 'N3D', 'IR', 'ASL', 'CNR', 'MP7', 'PBK', '4MU', 'RM4', '16U', 'EHD', 'RXD', 'G44', 'NLH', '701', '11R', 'PM9', 'OLP', 'Z8B', '3KC', 'TRF', '1F8', '3H1', 'OPT', 'AX7', '27U', 'M77', 'ISZ', 'TMS', 'PHB', 'OSL', 'L8P', '8JZ', '3HD', 'I1P', '175', 'SEE', 'AZ5', 'FHL', 'GVQ', '8DG', 'C7M', '895', 'LK1', '47C', '8NX', 'VAR', 'PIN', 'T48', 'GWJ', '3AL', 'H4B', 'DP9', 'LCE', 'AE1', 'JPZ', '2TY', 'CS3', 'AET', 'MX4', 'MHE', 'I2M', 'TDD', 'MND', 'HVA', 'LMQ', 'DHV', 'HTN', 'M2S', '359', 'FLL', 'MBQ', 'SC8', 'PBE', 'GVO', 'CSJ', 'COW', 'SGL', 'P2Q', '994', '2BD', 'R6A', '8TX', 'UMX', 'VSV', 'IAB', '6NA', '9PR', '1N9', 'FPR', 'BMP', '2A2', '9GP', 'CPA', 'OS', '5SE', 'BIM', 'AZK', 'RDD', '445', 'BSC', 'BTK', 'CLU', 'UC4', 'BGR', 'NRG', 'C4M', 'CTN', 'RCO', 'ZZL', '1SM', 'DGA', 'LAB', 'PXT', 'S2G', 'N5P', 'N5I', '26U', 'AD5', 'PRW', 'SKY', 'NBZ', 'X1E', 'ASO', '24X', 'G6Q', 'NIX', 'DUC', 'CNS', 'LY4', 'AP5', 'HDF', 'ZIT', 'NEG', '1N1', 'SBD', 'PAP', 'NYG', 'B17', 'XUL', '23V', 'DI6', 'TM2', 'HIU', 'N7P', 'DX4', 'OHT', 'RSQ', 'ISP', 'X37', '00I', 'G46', 'KSE', 'CH4', 'C2D', 'PC1', 'TAV', 'C1D', 'LEN', '0YG', 'QV4', 'PBZ', 'NCB', 'M09', 'L15', 'ZPR', 'TYC', 'TAO', 'VRV', 'MGT', 'ZYK', 'SCW', 'TRC', 'DNS', '770', '10U', 'C3X', 'M02', 'MEF', 'HXA', 'XTY', 'XCS', 'XGA', 'XAE', 'DMY', 'G1M', '276', 'IN2', 'SGV', 'S2M', '1GP', 'CAO', 'CIO', 'US3', '2AT', '458', 'GEN', 'BD2', 'DMX', '39A', 'MTZ', 'GVC', 'ADK', '2SC', 'ELY', 'SOY', '3HT', '2HT', '81A', 'JJK', '2LM', '5RM', 'LSO', '25D', 'TYO', 'LHL', 'PXA', '1LG', 'OLZ', 'CY0', 'DON', 'CR5', 'VX1', 'B2I', 'IBU', 'OBF', 'FE3', '128', '5AD', 'LL1', '170', 'OFF', '2DM', 'AN2', 'KSC', '517', '6CL', 'MOS', 'LG8', 'DQX', '741', 'IRB', 'URT', 'NCX', 'R96', '5FD', 'DUD', 'QRX', 'L02', 'IQP', 'TPC', 'AI3', 'ERY', '2PA', 'G93', 'RXB', 'TPZ', 'C2L', 'G2L', 'A2L', 'U2L', 'P34', 'VX3', 'TEL', 'M0H', 'DDD', 'SHF', 'I05', 'G2C', '4FW', '4F3', 'MPJ', 'TTT', 'B27', 'C19', 'NDS', '572', 'TA3', 'CIN', 'FMP', 'ZD6', 'F59', 'TIZ', '2JZ', 'PNT', 'XAN', '4PO', 'I04', 'CYJ', 'AAO', '23S', '56A', 'SLQ', 'FU1', '12D', 'DK5', 'HDE', 'X9Q', 'HMB', 'VAI', 'PQB', 'L0C', 'C34', '34P', 'AX3', 'G6A', 'GVK', '2AL', 'KN0', 'DHC', 'VGL', 'XJ1', 'P1P', 'GFH', '796', 'C99', 'HGL', 'AGD', 'LHC', 'CUD', 'KSF', '32S', 'MY3', 'DND', 'CJB', 'UD5', 'BG1', '4CP', 'T38', 'PRR', 'OXE', '1P6', 'CS4', '1IQ', 'URP', 'FC6', 'JMZ', 'CGP', 'BPH', 'LP6', '14A', 'N69', 'LW4', 'PAC', 'FMU', 'SC9', '5B1', '071', 'GX1', 'BID', 'L9N', 'M28', 'DKZ', 'PY3', 'JZZ', 'FK5', '17S', '144', 'FDG', 'GMS', 'FOM', 'ZRK', '5NC', 'YYA', 'ATZ', '26O', '2IA', 'Z2T', 'C95', 'ZYJ', 'DNO', '31U', 'GIR', 'CLY', 'UQ7', '92M', 'QLG', '23G', 'JNZ', 'ZRM', 'JJJ', 'ROF', '21N', 'GD', 'AMV', '39Z', '062', 'AM0', '2DP', 'PX4', 'THH', 'KN1', 'LY2', 'GYU', 'R55', 'KSM', 'BB2', 'I1N', 'B8L', 'PTL', '8I1', 'XIN', 'E1H', 'CWR', 'DXP', 'ARA', 'L9M', 'DZD', 'WIN', 'WFE', 'NSZ', 'DT5', 'C3Y', 'GFM', 'X2M', 'D3', 'XCT', 'XGU', 'XAD', 'XTH', 'CTR', 'FMY', 'B12', 'BTI', 'MFC', '5PB', 'G4D', 'FUD', 'NMT', 'I01', 'K6X', 'MDQ', 'KAN', 'YES', 'P29', 'DRN', 'P45', 'CXE', 'PIV', 'LPL', 'C92', 'KSH', 'NGT', '4LI', '34O', '152', 'OCB', 'D17', 'GAO', 'A5O', '24M', 'CTT', '3SC', 'HSQ', 'MLD', '240', 'GIL', 'C52', 'RRC', 'PYL', 'FC0', 'RGL', 'AAP', 'XPE', 'EFQ', 'RCY', '1TY', 'Z78', 'AZ1', 'SAS', '447', 'DX3', 'TH5', '3D1', '1CS', 'GL2', 'KPA', 'KCP', 'CLS', 'THN', 'REA', 'IDG', 'BDG', 'NEB', 'PCQ', 'AOG', 'GDU', 'IOM', 'LCP', '174', 'P3M', 'NEC', 'DP4', 'AP4', 'ITU', 'AZM', 'PKF', '645', 'LPT', 'PXG', 'ABM', 'MTE', '801', 'HG2', 'GAR', 'PLV', 'MXY', 'IDR', 'NMB', 'GTS', 'HPR', 'PNE', 'PCR', 'PGH', 'AON', 'CR3', 'SRY', 'MBD', 'HPT', 'NAA', 'AMI', 'PPG', '2BL', 'BEC', 'RPP', 'B2G', 'DHE', 'TMU', 'PAN', '10A', '334', 'MGM', 'BZI', 'LRB', 'DRB', 'JE2', 'AR', 'NTP', 'CK6', 'SHM', 'IPM', '711', 'REO', 'DSC', '0QB', 'N2P', 'OI1', 'I3P', 'FEM', 'PLO', 'TPV', 'LG2', 'SEI', 'ETP', 'DMJ', 'PGD', '6WO', 'F09', '25T', 'RTL', 'PID', 'CXF', 'HCA', 'CFN', 'CLF', '607', 'COT', 'PYE', 'XYD', 'BRH', 'AIC', '600', 'IOE', 'MAP', 'PFB', '907', 'KPC', 'MBO', 'DUT', 'PC', 'TCO', 'CK2', 'NAY', 'BR3', 'ANA', 'EQP', 'XCC', '9HX', 'IMR', 'KMP', 'VSO', 'ICL', 'NGP', 'UQ6', 'GTX', 'SDS', 'CMX', 'EQU', 'MTP', '787', 'CLI', 'DIZ', 'D6P', 'ACI', 'GLD', 'BLB', '429', 'MSC', 'ADX', 'ASV', 'BER', 'W35', 'HGB', 'CR6', 'CO4', 'GP8', 'DDF', 'SB', 'SBO', 'CUN', 'MCN', 'NG6', '2PH', 'CB4', 'GA', 'NPB', '4CO', 'RFP', '2MO', 'FTL', 'FUF', 'GPI', 'CEO', 'LI1', 'SQU', 'DNH', 'BBB', 'VXA', 'OMP', 'PDN', 'AIK', 'TME', 'L27', 'GDL', 'PHJ', 'BEP', 'DHT', 'FHP', 'NCZ', 'SIF', 'TBA', 'EAA', 'DEG', 'AMQ', '137', 'CH2', 'GUM', 'BX3', 'TPF', 'CK1', 'BDS', 'CH1', 'DM2', '6MP', 'LTL', '12P', 'MP1', 'YZ9', 'BPA', 'NHS', 'REP', 'IME', 'CBH', 'PHT', 'PPC', 'AXT', 'D12', 'OXD', 'HEG', 'ANH', 'AIP', 'PPZ', 'KAH', 'E4P', '2OS', 'BNI', '3PY', '2PG', '239', 'MID', 'U2G', 'L3P', 'L2P', 'GH3', 'PCW', 'TOY', 'TZ4', 'GL7', 'T44', 'XMI', 'GEL', 'FKP', '117', 'SPA', 'FMI', '0EM', 'PIY', 'PU9', 'I2P', 'IMH', 'DA7', 'RNP', '0IW', 'SFP', 'KEU', 'TSD', 'GTB', 'CDH', 'SAF', 'EMM', 'ADA', 'AL8', '3GP', 'FNG', 'LNQ', 'DK1', 'SPH', '0F7', 'SES', 'PU7', 'U0E', 'PTD', 'ESM', '772', '0GM', 'TCH', 'GNH', 'CLM', 'XDL', 'IGP', 'KDG', 'CHL', 'I06', 'CK5', 'NCN', 'RIL', 'MRC', 'DDH', 'FLU', 'IND', 'HH2', 'DPS', '485', 'TNE', 'CGF', 'AD3', 'ERT', 'BZF', 'Z5A', 'HCC', '4CA', 'LX1', 'TYK', 'HGX', 'TCT', 'T33', 'PPR', 'CVB', 'L24', 'CRA', 'CRN', '7HP', 'NLG', 'M7G', '433', 'PCT', 'OSU', 'MHC', 'PU5', 'ATR', 'OMX', 'ERE', 'MMQ', 'BBA', 'BCA', 'HHR', 'I4B', 'R36', '6FA', 'ABL', 'DCB', '821', 'AVG', 'LS1', 'IPD', 'CLR', 'DTN', 'AAD', 'BZP', 'MOB', 'IAD', 'NPL', 'PHP', 'P17', 'DEQ', 'KHA', 'E6C', 'MBG', 'TCL', 'ALR', '122', 'NT1', 'NT2', 'NA1', 'ADJ', 'TI3', 'GTM', 'ID5', 'HDA', 'ADB', '5AP', 'BLG', 'IYT', '154', 'NLA', 'KSG', 'PYQ', '397', '7RP', 'PDA', '444', 'XN2', 'PDE', 'GPB', 'KAI', 'D7P', 'FBT', 'RTB', 'IMF', 'NMQ', 'PVB', 'CEF', '3PP', 'CHM', 'GM2', 'TWT', 'PNN', '5RP', 'COH', 'EIC', 'CRT', 'PIM', 'GL9', 'FOK', 'PT1', 'REX', '2Z3', 'CLQ', 'CPU', 'DAU', 'SRL', 'PMS', 'BB1', 'PU3', 'PD', 'MUR', 'OHH', '843', 'PGX', 'PH3', 'NGO', 'PNM', 'ILO', 'CAD', '0DS', '142', 'PP8', 'UAA', 'ORX', 'DAN', 'OX1', 'TQ5', 'NDO', 'ITM', 'TEO', 'AHG', 'AP1', 'W84', 'ROX', 'NAQ', 'AMT', 'BFU', 'ASC', 'FCT', 'P0H', '0Q4', '9DI', 'DBT', '3PE', 'PAO', 'HSP', 'CYX', '1PN', 'TGF', 'YBT', 'A76', 'LDO', 'UP6', 'TET', 'L1P', 'L4P', 'ETE', '7AD', 'PBC', '172', 'FBQ', 'CSN', 'IMO', 'LPM', '0Z9', 'YF3', 'P23', 'YF4', 'SFM', '0GQ', 'RPN', '964', 'HLT', '127', '2AP', 'AMX', '653', 'EP1', 'GHA', 'GCV', 'PAF', 'FR9', 'DG6', 'DID', 'CP3', 'MGP', '232', 'RPL', 'M91', 'IN7', 'IMU', 'JEN', 'WY2', 'MNC', 'NOV', '5SD', 'PAV', 'GU4', 'GU6', 'GU0', 'GU5', 'GU8', 'GU9', 'GU1', 'GU2', 'GU3', 'GCT', 'GAD', 'THA', 'R79', 'HA1', 'QUM', 'TIA', '169', 'DDO', 'CP4', 'RHQ', 'RSS', 'SKP', 'UVW', 'CDG', 'NFS', 'DQU', '2FP', '2FA', 'CZH', 'PRZ', 'ANC', '3OL', 'GEG', 'RAC', 'W11', '4AM', 'PRT', 'UFM', 'FMF', '209', 'X04', 'EDT', 'DM4', '700', 'HC1', 'NSN', 'OCA', 'NAJ', 'NFN', 'NOS', 'RDF', 'FD3', 'DI5', 'HPY', 'TAL', 'DES', 'VS1', 'DRM', 'CAX', '5AN', 'CFM', 'U2P', 'AO1', 'MHZ', 'SPN', 'JST', '989', 'LEE', 'FDM', 'BTQ', 'AP6', 'BHI', 'TPD', 'HEP', 'BDA', 'CBB', 'U89', 'IM1', 'DP1', 'C48', '7NI', '7I2', 'PMO', '2BH', '34A', 'LIH', 'W42', 'CRB', 'MHB', 'ZAP', 'KT5', 'RBF', '2BM', 'AI2', 'IAG', 'AMB', '32P', 'TRL', 'RAS', 'APV', 'FG1', 'MOX', '4MZ', 'BIR', '6CP', 'GAA', 'P1A', 'EPH', 'FXN', 'EB1', 'PPE', 'PFZ', '213', 'CMI', 'AXP', 'PHW', 'IPT', 'TFP', 'PGM', '1BP', 'Y3', 'PI', 'H5P', '6NI', 'FTA', 'FBS', 'INL', 'GPX', 'CE2', 'NBA', 'ALH', 'CPM', 'NBN', 'PYZ', '3MT', 'PGU', 'A3S', 'EOT', 'CXB', '2HA', 'THU', 'PTX', 'XDN', 'TX4', '4NP', 'MK1', 'E10', 'MGR', 'MOT', 'MPM', '164', 'IOL', '270', 'CRM', 'DFX', 'BRD', '853', 'ARH', 'BNS', '44B', 'DCD', '5IC', 'MDC', 'LAM', 'SSC', 'WIA', 'A26', 'A2P', 'BEI', 'MMP', 'TDE', 'CMU', 'ACV', 'CDV', '901', 'SE4', 'U55', 'SBS', 'TMA', 'TR1', 'BPY', 'FRI', 'DNN', 'SNI', 'TAB', '5PA', 'PLU', '3CH', 'ASD', 'N', 'PFL', 'NTU', 'EMC', 'EMT', '34C', 'C26', '115', 'MA3', 'P4P', 'FR3', '146', 'MA4', 'ETR', 'URO', 'I5P', 'DM5', 'DAE', 'SSG', 'CCB', 'BJH', 'ND4', 'HO', 'T5A', 'LG6', '3NI', 'LII', '2AM', 'MH2', 'BNF', 'NBF', 'DAF', 'VER', 'PRM', 'IQZ', 'MNY', 'T16', 'O16', 'RIT', '785', 'CED', 'MC2', 'E4D', 'LVA', 'DBQ', '0ZL', 'BAA', 'PDP', 'HP2', 'TZC', 'DHQ', 'DRU', 'LRU', 'A88', 'LIO', 'HUP', 'TAP', 'HMU', 'AL7', 'CAB', 'NDR', 'SNG', 'GTR', 'LIS', 'PMZ', 'FPC', 'CIP', 'IN9', '1ZK', 'SJ1', 'ZTW', 'SM1', 'MYS', 'GSB', 'INH', 'CHC', 'HQC', 'NDE', '204', 'PAB', 'HEV', 'AZE', 'TSB', '300', 'PMT', 'STI', 'BGD', 'SRS', 'ST5', 'PPK', '45P', 'MNS', 'U66', 'SKF', 'PCG', '0DY', 'CPX', 'NCC', 'HNI', '2FU', 'IMT', 'PYB', 'DIB', '471', 'RCA', '1BH', 'SCV', '780', 'AKT', 'CTA', 'PUA', 'STG', 'CNH', 'SPO', 'SBX', 'U01', 'UKP', 'DAT', '1PT', 'FPE', '478', '2DI', 'DZF', 'NAE', '1Z7', 'TPT', 'FCX', 'RED', 'HUB', 'GSW', '7RA', 'IFG', '826', 'DP3', 'OPB', 'PTT', 'TS1', 'TS0', 'PLH', 'ZST', 'EIP', 'LXP', '1MZ', 'CP8', '097', 'GEQ', 'NDA', 'MF3', 'IDM', 'HEQ', 'TA1', '256', '7IN', 'AFN', 'TXL', 'DDP', 'TMQ', '3AR', 'D10', 'CZI', 'TND', 'HGU', 'IM2', 'QUS', 'CTY', 'BLS', '0ZR', '4HA', 'FNP', 'PBM', 'BRC', '493', 'IAV', '017', 'NSC', 'IPP', 'MM3', 'HG9', '666', '2NH', 'APU', 'IET', 'YSH', 'GPJ', 'S3P', 'TQD', '6IN', 'BHA', 'VAA', 'SNR', 'PDX', 'PAL', 'ITR', 'A77', 'BGL', 'MQ8', 'EPB', 'PGS', '5PH', 'MCA', 'NAU', 'NAV', 'XIM', '3BT', 'NDC', 'CYH', 'CAG', 'FR1', 'XMF', 'H4P', 'PU2', 'OXN', 'MLN', 'IPB', 'TPX', 'SDC', 'HAR', 'DRF', 'OMD', 'FMC', 'BIP', 'NPG', 'W05', 'QNO', 'MQU', '9MG', 'FUP', 'NAX', '0I5', 'AST', 'NHB', 'CHT', 'FAF', 'GTG', 'FII', '12H', 'PNA', 'BBT', 'PU0', 'MN3', 'DHM', 'FLP', 'STC', 'FPY', 'TI2', '5OB', 'NTZ', 'AY1', '3NH', 'MAQ', 'PPD', 'FAM', 'GPP', 'H2B', 'EDE', 'DFV', 'EBS', 'BEM', 'MAV', 'LGU', 'IMG', 'BTS', 'GUP', '0QF', 'AHE', 'CMS', 'FOS', 'SM2', 'INE', 'IDP', 'EPA', 'TK4', 'AFP', 'PNP', '9AR', 'LPR', 'ENO', 'BFI', 'RP1', 'HC0', 'BEG', 'MCT', 'DAV', 'SND', '14W', '118', '5IN', 'TFM', 'FPH', '654', '679', 'NTA', 'CCE', 'KMT', 'EQI', '818', 'C15', 'CNL', 'ANN', 'AND', '750', '977', 'BEJ', 'PSR', '852', 'PPS', 'VTQ', 'SUM', '869', 'BDD', 'ECO', 'G6D', 'VO3', 'ACG', 'PON', 'DEF', 'HHA', 'ENH', 'DME', '340', 'ATS', '7DG', 'DGG', 'WW7', 'PQ0', 'I11', 'NPE', 'X7O', 'BGX', 'CL4', 'LMS', 'ETN', 'MSF', 'C78', 'IMA', '545', 'CFP', '961', 'OPS', 'FMA', 'SS1', 'MYA', 'MIM', 'P13', 'EAL', 'C2O', 'UDX', 'PNG', 'XIL', 'AOA', 'BV4', '656', 'GU7', 'RHP', 'HXC', 'FCB', 'F2B', 'FMB', 'TOP', 'D16', 'RSO', 'FRE', 'PCB', 'ARR', 'ZEN', '101', '0HG', 'AH1', 'PBG', 'PY6', 'DEX', 'FIP', '5MP', 'BFA', 'INS', 'NBC', '207', 'VS2', 'F6B', 'GYP', 'UCM', 'RDC', 'PIC', 'ELP', 'SF3', 'PNS', 'MYG', 'SB5', 'OTS', 'OTR', 'HEH', '106', 'T4A', 'CO8', 'ENC', '4BT', 'BP1', 'TRH', 'IMM', 'NAT', 'IPL', 'DIG', 'GDN', 'DOM', 'MBP', '0Z4', 'C2N', 'NAK', 'BDM', 'G2F', 'NFG', 'TYV', 'XLC', 'C5G', 'PCP', 'MC9', 'OPD', '4PT', 'CXL', 'RU7', '383', 'OIN', 'GP6', 'AML', '3AG', '806', 'PEJ', 'MAO', 'BHG', 'LNL', 'GRO', 'TBO', '2A6', 'FNE', 'TTO', 'FRG', 'DQO', '0PP', 'BZA', 'AIG', 'VS3', 'SOG', 'HAZ', 'PUZ', 'QUN', 'VDX', '3AI', 'KTN', 'PVH', 'PCL', '941', 'TOA', 'TOC', 'TOB', 'HQQ', 'DBF', '0ZQ', 'C1O', 'M2C', '183', 'BRP', 'DPD', 'PCN', 'BP2', '13S', '13R', '9OH', '11O', '6PG', 'P3G', 'ACQ', 'UP5', '2DG', 'FLV', 'TNT', 'GBS', 'ETF', 'ADL', '0FP', 'DBG', 'DTL', 'CR4', 'TDA', 'AAG', 'PHY', '0H8', 'SBA', 'ILA', 'NBE', 'MIC', 'DTO', 'FR7', 'GAG', 'MNX', 'ACD', 'FFO', 'DCM', 'PA7', 'W03', 'CNN', '161', 'INK', '778', 'MEV', 'BPG', 'FD4', '537', '1MC', 'CTD', 'NGH', 'HEE', 'TZD', 'GIO', 'G26', 'LYA', '084', 'EPX', 'MX1', 'CCO', 'GKR', 'PYN', 'OPA', 'OXQ', 'PML', 'ANB', 'AGN', 'IYG', 'IC1', 'DCU', 'SL1', 'MJI', '762', 'SAD', 'KT3', 'GLO', 'DHY', 'GUA', 'PTE', 'LUM', 'LCO', 'ILT', 'H2A', 'L79', '570', '185', 'CBP', 'DUX', '802', '0DO', '2MP', 'INP', 'GTK', 'SD8', 'SPF', 'XV6', '074', 'PDM', 'VDY', 'KR', 'MF2', 'DI2', 'MSD', 'MNO', 'CDC', '292', 'CN2', 'FAE', '33P', 'UAP', 'HYC', '675', 'FNS', 'NBU', 'IH5', 'AR1', 'R01', 'DAX', 'CLL', 'E12', 'BIK', 'C20', 'BC', '113', 'SRM', 'BM2', 'RAI', 'DI4', 'GAP', 'UNA', 'ST4', 'BBH', 'PPQ', 'PAD', 'MYC', 'DR1', '4CM', 'L04', 'NLX', 'FSF', 'BN1', 'DRG', 'NBD', 'BG5', 'MCR', 'RWJ', 'BZB', 'GUD', 'HF5', 'HF3', 'PHF', '0PO', 'SHR', '2CM', 'NAF', 'HMN', '1PM', '4AP', '0DB', 'KBG', '2Y3', 'HBU', 'BGN', 'TMC', 'MPH', 'BG3', '2NP', '3PC', 'ST6', 'MAF', 'MBF', 'NIN', 'CCH', 'DAK', 'IN0', 'DIA', 'NDT', 'IK2', 'A79', 'M6T', 'PRA', '2CP', 'DIR', '1GL', 'ARI', 'CPH', 'CDR', 'ERI', 'FBI', '667', 'SSE', 'MHN', '4NC', 'SN2', '108', 'TIN', 'IRP', 'PMM', 'RAP', 'DMD', 'SGP', 'ONM', 'URF', 'A8B', 'PSB', 'HCO', 'LY1', 'HPS', 'A45', 'FEL', 'MPN', '141', 'IPZ', 'PRC', 'KIF', 'FCR', 'GCD', 'LVS', '214', 'MSP', 'ETI', 'NTH', 'IBM', 'DO2', 'UN4', 'AQ4', 'Q82', 'LHY', 'LAK', 'DBM', 'TDS', 'OPC', 'MFU', 'ATF', 'CPO', 'FQP', 'PZQ', 'CBD', 'DQN', 'TI1', '234', 'BHH', '3EP', 'MXA', 'BMU', 'HAX', '1PU', 'GTZ', 'TAQ', 'D1B', 'MCI', 'OCP', '4MV', 'TPE', 'ATY', 'D13', 'LS4', 'R17', 'BAK', 'AL3', 'IP2', 'FFB', 'A32', 'UI3', 'B96', 'R18', 'D4T', 'D4D', 'DX9', 'SIM', 'DP2', '2MZ', 'ODP', 'CEP', 'HTL', 'SB4', 'ESR', 'DDT', 'FEN', 'DIC', '2FD', 'R88', 'CPF', 'HSD', 'OOA', 'N4B', 'STL', 'COY', 'COB', 'CBO', '2PB', 'CR1', 'MOC', 'HGP', 'NEQ', 'CLP', 'UMG', 'IU5', 'S2C', '1Z1', 'MUP', 'MIL', 'KEG', '681', 'PSI', 'AGP', 'SON', 'MAY', 'INI', 'KTP', 'ALL', 'HUX', 'PAM', 'TN4', 'CIL', 'MLP', 'SHI', 'MIT', '588', 'ACX', '24T', 'ROI', 'CPR', 'BMQ', 'DYB', 'IN3', 'CIU', 'TMM', 'ESA', 'DAC', 'BAY', 'MTH', 'BIH', '655', 'OAP', 'DDI', 'DHR', 'RAE', 'C1P', 'MOP', 'HBO', 'FSP', '2Y2', 'MTI', 'RSA', 'CEL', 'PUR', 'PUX', 'IBS', '497', 'PYD', 'P25', 'MGC', 'BV2', 'J77', 'ZNH', 'KHP', 'GNT', 'WCC', 'UI1', 'TDX', '094', 'HAM', 'FLF', 'ISA', 'SAE', 'MBR', 'CBT', 'PY1', 'DER', 'PHH', 'NGV', 'PBP', 'DDY', 'IS2', 'BAU', 'WRA', 'CEQ', 'ZIP', 'ADS', 'THS', 'TPM', 'VWW', '936', '2BF', '114', 'NCT', 'N25', 'HCG', 'STB', 'CRI', 'CPW', '0EF', 'GPS', 'NPP', 'IVV', 'GDR', 'CP0', 'BM6', 'PNU', 'CK8', 'IPO', 'AL4', 'PBN', 'AAS', 'NTD', 'SUW', 'NPR', '846', 'BM9', 'PMD', 'SPC', 'LIP', 'CNB', 'CN1', 'CNF', 'CIA', 'AN1', '87Y', 'NON', 'TNC', 'UD2', 'PET', 'HE3', 'M1A', 'KSA', 'PAX', 'PR2', 'E2P', 'BP6', 'PDT', 'BCZ', 'FID', 'AS1', '336', 'NA4', '15B', 'J15', 'AMR', 'BZM', '794', 'SPL', 'GND', 'PEH', 'ALJ', '2ZS', 'J12', 'ASE', 'MTD', 'TE', 'HDN', 'KHO', 'MNT', '4BZ', 'CFO', 'IFL', 'R99', '3AA', 'BZS', 'BMD', 'KAM', 'CGQ', 'ISC', 'DFB', '104', 'HTP', '968', 'XMH', 'PHV', '140', 'R03', 'IUR', 'EU3', '1C5', 'FSB', '9DA', 'CRZ', 'BPO', 'CTB', 'CBA', 'CZM', '145', 'PSN', 'DVC', 'BML', 'MZM', 'MD2', 'ZEM', '3FM', 'PMY', 'PSO', 'HEL', 'COG', 'VD1', 'I52', '2FL', 'R71', 'F23', 'GAS', 'EG3', 'RR6', 'SM3', 'ATQ', 'MSI', 'ST1', 'AQO', 'IOC', 'IOT', 'GLG', 'ARM', 'ARV', '1ZG', 'AE2', 'DOB', 'PA2', 'P24', '42B', '182', 'MCG', 'ING', 'COP', 'HBY', 'SLT', 'TX5', 'HAS', '0QN', 'TOE', 'PAJ', 'C2C', 'QSI', 'PYF', 'TAA', 'FBU', 'APS', 'ARP', 'MLM', 'ABI', '1UN', 'BR5', 'RCL', 'IDI', 'HST', 'RIS', 'EMD', 'TBN', 'DTI', 'BA3', 'UGA', 'AXL', 'BPJ', 'BRB', 'ONP', 'MNP', 'DNC', 'COS', 'PSF', 'POB', 'LKA', 'NPI', 'SCO', 'BHS', 'DIH', 'NF', 'D24', 'NMP', '5MD', 'FRC', 'NFL', '418', 'PRH', 'PU4', '2Z0', 'ABD', '24B', 'HYG', 'ICP', 'TZP', '171', 'PEQ', 'A70', 'LGZ', 'MGU', 'SHH', 'OLO', 'FTH', 'AMH', '4PN', 'GB', 'DST', 'T29', '155', 'IN', 'RXP', '335', 'VSC', 'DHS', 'CK3', 'BZR', 'DMB', 'PMC', 'C2P', 'ICU', 'DM7', '1RB', 'BUQ', 'OAI', 'INQ', '442', 'UCL', 'DFL', 'DIT', 'XK2', 'DNQ', 'LS3', 'DHH', 'PB1', 'TPA', '272', 'MAB', '1CP', 'NVI', 'PLS', 'CDN', 'AAA', 'AOP', 'PNQ', 'MTT', 'P2Y', 'PRX', 'EPG', 'AAH', 'NFH', '160', 'TZL', 'MOD', 'ENQ', 'BZQ', '0PQ', 'W71', 'R37', 'FWD', 'J80', 'SKM', 'FA3', 'RBZ', 'BPT', 'YMA', 'CE', 'IXM', 'DAQ', 'SOA', 'CNO', '847', 'RBU', 'DEB', 'CAZ', 'SBT', 'PE2', 'PG2', '3DG', 'SB1', 'CFX', 'OUT', '0ZY', '773', 'SN1', 'INV', 'HFM', 'IHG', '2Z4', 'P1C', '950', 'CLW', 'BTT', 'ADV', '0ZJ', 'MUG', 'RPD', 'TPB', 'AP3', 'EP2', 'AHX', 'XYH', 'T1P', 'T2P', 'T5P', 'T4P', 'IPY', '135', 'HFL', 'IR3', 'MAH', 'AT1', 'MLR', 'GTA', '0D3', 'NHP', 'T10', 'NMD', 'PCM', 'TCK', 'AIQ', 'LAF', 'OQB', 'OCV', 'GLE', 'UFG', 'PHN', 'PIW', 'FAS', 'HWG', 'R46', 'SSD', 'GAT', 'TEA', 'CO1', 'FRQ', 'URS', 'MUS', 'HEB', '3ID', 'RDR', 'A24', 'N1T', 'FON', 'DMI', 'INA', '629', 'PMB', 'PA3', 'H3S', 'BRW', 'AOE', 'TNF', 'I10', 'DEN', 'BLI', 'MOF', 'UND', 'TRD', 'XDP', '4HM', 'BH0', 'SU9', '815', 'FTI', 'DVR', 'BRE', '5AS', '2SA', 'IMY', 'PYG', 'CBN', 'APX', 'DMV', 'XMA', 'GMC', 'GR3', '1PB', 'OXP', 'BRL', '434', 'NOG', 'TPY', 'SC4', 'CMA', 'N2M', 'UHD', 'CRD', 'DNT', '3AN', 'PCY', 'V', 'IMB', 'MER', '134', 'VIG', 'G37', 'MS3', 'NDM', 'PL1', 'CR', 'SSO', 'WAY', 'ADH', 'HI5', 'SMS', 'TMH', 'GZZ', 'EPT', 'GTD', 'YSA', 'TBR', 'NPF', 'A15', '1NH', 'NES', 'DHZ', '0HH', 'ANQ', 'APQ', 'IWD', 'P6C', '2AC', 'DPB', 'CYE', 'ISN', '460', 'PYS', '105', 'TRZ', 'MAW', 'AHM', 'YPA', 'MAT', 'DPA', 'DAG', 'TMR', 'DLF', 'MTC', '0A0', 'EEB', 'TRR', 'AFI', 'C4P', 'RTC', 'DPT', '1QQ', 'R4K', 'AKK', '688', 'FLM', 'DPY', 'DRP', '0ZZ', 'CES', 'KEF', 'CYZ', 'XIF', 'A12', 'IN4', 'RTA', '783', 'HGI', 'RTR', 'HLE', 'RIN', 'MBB', 'HDM', 'GAN', 'FHC', '262', 'TBP', 'DED', 'BRZ', 'EGT', 'LPE', 'TMT', 'GLR', 'DTM', 'NHO', '0P0', 'CP5', 'VIO', 'KMB', '0E7', '197', 'PAH', 'E20', 'MM4', '3AD', 'NDH', 'MYD', 'CXU', 'HAB', 'TCZ', 'WY4', 'VA3', 'V7O', 'NAB', 'MOG', '0Z0', 'SOX', 'RRS', 'HPL', 'OGG', 'PSQ', 'OVA', '2PU', 'FRA', 'MO7', 'OMO', 'BMS', 'MTM', 'TST', 'BDK', 'NMH', 'AL5', '156', 'GNB', 'IP1', 'AL2', 'POD', 'E09', 'IL0', 'NLP', 'R64', 'SWF', 'TDR', '4AB', 'THT', 'YBH', 'MTS', 'CEI', 'TEI', 'CCV', 'AHI', 'AM1', '0ZX', '158', 'SH4', 'ZAR', 'PPJ', 'GP3', 'VG1', '972', 'HOP', 'D1L', '5OP', 'SBR', 'BCG', 'DDC', 'PRY', '3GR', '151', '5UD', 'RAO', 'GBD', 'EBW', 'GE1', 'GE2', 'GE3', '061', '0EG', 'BHP', 'UPA', '2FH', 'NHD', 'PSG', 'ITP', 'EBG', 'MPC', 'WAC', 'DEO', 'TCC', 'NPY', 'FAA', 'DQA', '109', 'MYY', 'BOM', 'BIS', 'T80', 'A85', 'W01', 'LMZ', 'BDF', 'UDH', 'LDT', 'DQH', 'DH2', 'S70', 'EPY', 'ME2', 'DXE', 'AEN', 'LPG', 'FD1', 'INW', 'BRR', 'BL5', 'PY4', 'TQ4', 'APZ', 'SMG', 'BED', 'M2M', 'P5A', 'BCH', 'PRI', '791', '0ZG', 'PCX', 'D34', 'PG3', 'TCN', 'MSQ', 'A8N', 'M1C', '693', 'ESP', 'DCO', 'UR2', 'AL1', 'AG7', 'RLP', 'PHS', '254', 'ENP', 'L34', 'BEK', 'DSA', 'PHG', 'ET', 'DCA', 'SCT', 'SPV', 'NMX', 'ST2', 'GTY', 'BNR', 'NXN', '0HZ', 'NDD', 'PH1', 'T6P', 'HPF', 'DOS', 'LOS', 'J78', '892', 'FUX', 'OST', 'INM', 'LTN', 'EDR', 'WRS', '0Z3', 'CA5', 'NAR', 'AGB', 'PCO', 'DXG', 'BES', 'BDI', 'DTC', 'OKA', 'CYY', '544', '2IP', '299', 'XLD', 'BD1', '564', 'IDX', 'AMC', 'AOM', 'C3S', 'DPC', 'SHV', '107', 'BE1', 'UN2', 'PMA', 'NH1', 'SIL', 'RMN', 'SMN', 'CTO', '486', 'GC7', 'TBT', 'FLA', '5CA', 'PBO', 'SDZ', 'MXG', 'GUL', 'FBA', 'NOP', 'TQ6', 'R56', 'SHU', 'ZEA', 'GIP', '0Z1', 'MYP', 'CDT', 'LYB', 'UI2', 'TQT', 'BHB', '2MC', 'ESO', '1DM', 'SMO', 'GS1', 'BDB', 'GCW', 'BV1', 'SUA', 'XMD', 'I58', '974', 'EG2', 'CER', 'FAP', 'FCP', 'BYS', '103', 'MLG', 'NSP', '2AF', 'MGS', 'GDA', 'FLG', 'LAZ', 'BH7', '216', 'AV2', 'SI1', '5BN', 'CA3', 'PEM', 'FEP', 'LFA', 'PAI', 'CFC', 'HOM', 'BEL', 'MFB', 'ACC', 'HNA', 'EG1', 'CS8', 'INF', 'FKA', 'SDK', '178', 'GAF', '147', 'BHC', 'BP3', 'AYM', 'GW5', 'JAN', '0ZW', '623', 'CHO', 'PBR', 'SNA', 'PR1', 'LIG', 'AB1', 'IHE', '8IN', 'SH1', 'FRJ', 'SQA', 'TDB', 'CEN', 'PRB', 'TMZ', 'S58', 'CBF', 'P2C', 'DMZ', 'AZ2', 'CPQ', '001', 'EJT', 'FOH', '3HA', 'SBP', 'NFE', 'CAI', 'SPR', 'TDT', 'BLT', 'NID', 'ALI', 'PPB', 'RFL', 'OSS', 'MYE', 'BRS', 'MGB', 'N3T', 'FR0', 'CDX', 'BTH', '153', 'B9A', 'RE9', 'EP0', 'CVI', '44D', 'SG2', 'TFK', 'TAX', '312', 'XN1', 'OBA', 'BND', 'R19', 'LY3', 'ATU', 'AZR', '5MB', 'D76', 'ISQ', 'GPC', 'B1V', '878', 'L08', 'UVC', '01S', 'RUT', 'DSS', '739', 'SP1', 'ROP', '0ED', 'FEA', 'OFL', 'DDA', 'DDL', 'MDA', 'DXB', 'MQD', '3CY', 'PS5', 'S27', 'SSB', 'BNL', 'DSO', 'BV3', 'EP', 'PM1', 'FD2', 'KPL', 'HE1', 'AO5', 'XN3', 'AEP', 'TSN', '2BN', 'TBE', 'MDD', 'ETS', 'RHC', 'PP7', 'ZAM', 'FRZ', 'JY', 'ZK5', 'HBC', '4HC', 'GL4', '0P1', 'BSB', 'OBE', 'RA2', 'IDA', 'NHM', 'SB6', 'MDR', 'KNC', 'PF3', 'PFE', 'MPL', '8BR', 'BUL', 'IN6', 'P90', 'GCH', 'HTQ', 'AKR', 'TAZ', 'PU8', 'BSI', 'MYT', 'IOP', '616', 'BEH', 'B3N', 'AR3', 'DEC', 'APG', 'PU1', 'FR2', 'S80', 'PRD', 'CEH', '186', '5FA', 'BBM', 'IDN', 'NMA', 'NM1', 'PTU', 'A78', 'APT', 'ZEB', 'HYF', 'BG4', '822', 'TBI', 'CBQ', 'ZMR', 'GNR', 'M3C', 'THD', 'W37', 'AA2', 'DBR', '35G', 'NCH', 'AZC', 'BAX', '5EA', 'L75', 'AMN', 'MTB', 'CGS', '876', 'KWT', 'AC2', 'GCR', 'FMD', '123', 'GM6', 'HPI', '0PX', 'BNE', 'ABP', 'DFN', 'BC1', 'SPB', 'HWD', 'DEM', 'XME', 'WRB', '303', 'C60', 'FU2', 'K57', 'GMN', 'ANO', 'I40', '2PC', 'Q50', 'RE1', 'FRH', '2EZ', 'TBC', 'ISF', 'MHA', 'FTP', '116', 'ESX', 'NTS', 'FA6', 'EH5', 'TTA', 'CAH', 'BFL', '0ZN', '761', '9AP', '184', 'RMA', 'BRF', 'NHR', 'DKT', 'L37', 'NTM', 'BOA', 'OPP', 'BPM', 'NZQ', 'DOG', 'DUO', 'MFP', 'MFQ', 'OTA', 'PCV', 'FYA', 'HH0', 'NFT', '687', 'CZZ', 'UNI', 'AM4', 'GTL', '515', 'ICA', 'MEI', 'CTS', '678', 'ATC', 'DGX', 'SUB', 'FMS', 'FA5', '133', 'COC', 'MRK', 'FDQ', 'CHQ', 'CYU', '998', 'RUN', 'ABZ', 'FM2', 'FM1', 'RAB', '159', 'INN', 'BSP', '426', 'KNI', '0QH', '100', 'GGD', 'I48', 'DDN', 'TNR', 'TBH', 'KOS', 'A3B', 'NIP', 'GOM', '4TP', 'DCW', '4BR', 'NA3', '03D', 'MTQ', 'MON', 'LS2', 'TOT', 'COQ', '181', 'LGC', 'SYM', 'ABB', '4PA', 'GAC', '0E8', 'FEX', 'G24', 'CEM', 'TEM', 'OIR', 'NGZ', 'DBS', '991', 'BFS', 'PU6', 'IPU', 'DBD', 'NYP', 'NCR', 'SU2', 'ID2', 'SUP', '0A8', 'X2F', 'HPN', 'VS4', 'C2G', 'LAP', '149', '0QI', '9AM', 'CUM', 'KET', 'HEN', 'CP9', 'FL2', 'FFC', 'DCH', 'CUO', 'SPW', 'LYD', '2CA', '4TA', 'CCI', 'P10', '703', 'NXA', 'FEC', 'FCN', 'R02', '880', 'NOE', 'XMJ', 'CMB', '4FA', 'AZD', '803', 'BEB', 'EPD', 'EMA', '2LP', 'IOB', 'LPC', 'CXA', 'PN1', 'AEM', 'APD', 'DFA', 'H37', 'MDE', 'WSK', 'RTP', 'SHY', 'XMK', 'HMO', 'DRC', 'MNA', 'MPA', 'DAJ', 'PRE', '745', 'CIB', 'BU2', 'BB3', '8IG', 'FUG', 'MPE', 'BJP', 'DUP', 'HCP', 'HMD', 'IMZ', 'MBV', 'AZN', 'NNH', 'CLZ', 'D35', 'RS7', 'NCG', '1AN', 'RMB', 'GL5', '797', 'ST3', 'NIR', 'AL6', 'GET', 'DLG', 'OXS', 'DP7', '0EO', 'NOC', 'INX', 'CLC', 'MS1', 'PZN', 'SPQ', 'GP2', 'PIO', 'RNT', '49A', 'GCG', 'IQA', 'PC5', 'PE6', 'EMU', 'P2P', 'PIS', '0E6', 'PAT', 'AZG', '65B', '0IU', 'HEY', 'PIQ', 'VD2', 'PFS', 'SAN', 'SBB', 'Z34', 'MZP', 'ADY', 'I17', 'ALB', 'B4G', 'IDH', 'QNC', 'DIQ', 'NTB', 'BDC', 'EMB', 'MEC', 'HDC', 'LUT', 'XAT', 'NEX', 'HIF', 'NPC', 'LYL', 'FYX', 'SCM', 'TQ3', 'DNF', 'CUZ', 'IDU', 'BTD', 'PFC', '111', 'TZ5', 'IOF', 'AAN', 'BWD', 'RU', 'AHF', 'APL', 'APF', 'DOT', 'MYX', 'PLX', '2ZF', '2GL', '1FH', 'AKA', 'SUD', 'AP', 'WO5', 'GW3', '9PP', '4CB', 'TLM', 'IBA', 'PHK', '0ZB', 'ABT', 'PUY', 'RNS', 'SG1', 'FXY', 'IHJ', 'WRR', 'BMV', 'C11', 'JEF', 'RPA', 'TEP', '2C5', 'NYL', '2Y4', 'P27', 'SBN', '338', 'PE7', 'IBG', 'PEZ', '903', 'LVG', 'BW2', 'LKS', '4QB', 'IVS', 'R23', 'XBP', 'IBB', 'ATA', 'COL', 'DRI', 'DXA', 'CPZ', 'LOR', 'ACZ', 'INZ', 'GMY', 'GMM', 'PSD', 'IDD', 'CRH', 'CG', 'ZYA', 'TYP', '0D6', 'MDN', '0FG', 'RR1', 'THK', 'FTB', '358', '124', 'AI1', 'ALP', 'LP1', 'CLK', 'LG1', '139', 'FL9', 'RDI', 'VX', 'MBS', '394', 'IOA', 'TZA', 'NOM', 'TNL', '858', 'B7G', '2IN', 'SMB', 'ARL', 'PYM', 'CDU', 'ZFB', 'HES', 'AAM', 'NPA', 'DEA', 'AO2', 'A3M', '110', 'HYB', '7A8', 'FPI', 'HAV', '3MP', 'MSB', 'TTX', 'ABH', 'V36', 'DI', 'ACS', 'B3I', 'SWA', '3HC', '680', 'FR5', 'E3G', 'TCB', 'FOP', 'DBI', 'DSI', '997', 'OAL', 'LYC', 'DET', 'SPI', 'BCP', 'TSA', 'BZU', 'TSX', 'RRR', '326', 'VIT', 'PND', 'BAF', 'OFO', 'UZ9', 'GA2', 'PP9', 'SU1', 'RRP', 'PTI', 'RPR', '9TA', 'FL8', 'EMO', 'XMC', 'SIH', 'CNE', 'BPD', 'RFA', 'RFB', 'GBC', 'AJ3', 'MQ1', 'SK1', 'K21', 'FA1', '9HP', 'BLN', 'STN', 'BUM', '1DA', 'DTX', 'CND', 'BIA', '2FM', 'IKT', 'N5B', 'IMI', 'FR8', 'USQ', 'TPI', '138', 'CUB', 'RAD', 'MNQ', 'ERG', '150', 'FRB', '129', 'I3N', '2HC', 'IFP', 'MDM', 'SBZ', 'P2S', 'Q2Y', 'POS', 'DO3', 'XMB', '587', '2NO', 'CMM', 'INJ', 'PPO', 'RQ3', 'P14', 'BRT', '3AP', 'RIF', '5NI', 'PY5', 'RH1', 'SAB', 'TFI', 'W02', 'CZN', 'HBN', 'AKP', 'AHU', 'KH1', '219', 'DHP', 'RDL', 'U49', 'IHI', 'FDI', 'FAL', 'FBL', 'TSU', 'IB2', 'ASR', 'MG7', 'AAT', 'SD2', '354', 'FXV', 'PZO', 'BVD', 'SLF', 'P12', 'R04', 'VK3', 'SP2', 'BR4', 'LXC', 'HE6', 'BLV', 'U3H', 'DI3', 'U2F', '3IN', 'NPN', 'HBS', 'HBR', '984', 'CON', 'FRM', 'BZC', 'TC4', '760', 'CHR', 'PTS', '580', 'TTE', 'TBD', '1NB', 'E97', 'G23', 'FOC', 'LS5', 'HAG', 'ART', 'BUB', 'NEO', 'ILP', 'LOX', 'ZES', 'BYP', 'BWP', 'HF1', 'INB', 'LCS', 'A4P', 'GBP', 'DAA', 'AUR', 'AYD', 'RL2', 'AU3', 'PXM', '4HP', 'IBP', 'AL9', 'COJ', '136', 'HTA', 'C3P', 'RSD', 'XMG', 'CP6', 'CDY', 'BVP', 'VDZ', 'M90', 'HBD', '5YL', 'DPU', 'SB3', 'BDU', 'MA2', 'MA1', 'APO', 'P11', 'S57', 'AF1', 'ALD', 'BBZ', '4PB', 'CC1', '965', 'TS5', 'PFA', 'CC0', 'CXT', 'LDM', 'REY', 'CB1', 'KOJ', 'RBE', 'FR4', 'BEE', 'ROM', 'NBP', 'MR2', 'LAX', 'HDI', 'LAD', 'ELD', 'DPZ', 'GEO', 'PII', 'MPI', 'FUR', 'SBE']

#######################   ClashRecognizer   ##########################
#from Moderna Constans
ATOM_RADII = {'C' : 1.72, 'N' : 1.55, 'P' : 1.8, 'O' : 1.52, 'S' : 1.8, 'H': 1.2, 'I': 1.98}
# All radii above were taken from 
# http://geometry.molmovdb.org/files/
# libproteingeometry/data/NucProt.atom-defs.dat

#SEARCH_RADIUS = 2.1 # max sum of two radii


################### Allowed chain names ############################


ALLOWED_CHAIN_NAMES = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","W","X","Y","Z",\
                       "a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","w","x","y","z",\
                       "1","2","3","4","5","6","7","8","9","0","-","+","_","=","~","`","!","@","#","$","%","^","&","*","(",\
                       ")","{","}","[","]","|"]

################################################################

####################  Molecular Weights ################################

## For Center of Mass Calculation.
## Taken from http://www.chem.qmul.ac.uk/iupac/AtWt/
MOLWEIGHTS = {
             '?' : 0.0,     'H' :  1.00794,   'C' : 12.0107, 'N' : 14.0067,
             'O' : 15.9994, 'P' : 30.973761,  'S' : 32.065, 'I': 126.9}

##############################################################################


################  Volume Calculation   ######################################
AA_WEIGHTS = {"A":71.08,  "C":103.14, "D":115.09,
              "E":129.12 ,"F":147.18, "G":57.06,  "H":137.15,
              "I":113.17 ,"K":128.18, "L":113.17, "M":131.21,
              "N":114.11 ,"P":97.12,  "Q":128.41, "R":156.20, "S":87.08 ,
              "T":101.11, "V":99.14,  "W":186.21, "Y":163.18, "X": 110, "U": 110}
##The molecular weights above are those of the free acid and not the residue ,
##which is used in the claculations performed by the Peptide Properties Calculator.
##http://www.basic.northwestern.edu/biotools/proteincalc.html
##Subtracting an the weight of a mole of water (18g/mol) yields the molecular
##weight of the residue. The weights used for Glx and Asx are averages.


#@TODO: weights for X (modified) residues are mean value from known mol weight for common nucleotides
DNA_WEIGHTS = {"A": 313.21, "C": 289.18, "G": 329.21, "T": 304.2, "X": 308.95, "U": 306.17}

RNA_WEIGHTS = {"A": 329.21, "U": 306.17, "C": 305.18, "G": 345.21, "X": 321.44}

############################################################################

############## Residue names (three letter to one letter conversion) ##############

#from Bio.Data.ScopeData in BioPython

to_one_letter_code = {
    '00C':'C','01W':'X','02K':'A','03Y':'C','07O':'C',
    '08P':'C','0A0':'D','0A1':'Y','0A2':'K','0A8':'C',
    '0AA':'V','0AB':'V','0AC':'G','0AD':'G','0AF':'W',
    '0AG':'L','0AH':'S','0AK':'D','0AM':'A','0AP':'C',
    '0AU':'U','0AV':'A','0AZ':'P','0BN':'F','0C ':'C',
    '0CS':'A','0DC':'C','0DG':'G','0DT':'T','0FL':'A',
    '0G ':'G','0NC':'A','0SP':'A','0U ':'U','0YG':'YG',
    '10C':'C','125':'U','126':'U','127':'U','128':'N',
    '12A':'A','143':'C','175':'ASG','193':'X','1AP':'A',
    '1MA':'A','1MG':'G','1PA':'F','1PI':'A','1PR':'N',
    '1SC':'C','1TQ':'W','1TY':'Y','1X6':'S','200':'F',
    '23F':'F','23S':'X','26B':'T','2AD':'X','2AG':'A',
    '2AO':'X','2AR':'A','2AS':'X','2AT':'T','2AU':'U',
    '2BD':'I','2BT':'T','2BU':'A','2CO':'C','2DA':'A',
    '2DF':'N','2DM':'N','2DO':'X','2DT':'T','2EG':'G',
    '2FE':'N','2FI':'N','2FM':'M','2GT':'T','2HF':'H',
    '2LU':'L','2MA':'A','2MG':'G','2ML':'L','2MR':'R',
    '2MT':'P','2MU':'U','2NT':'T','2OM':'U','2OT':'T',
    '2PI':'X','2PR':'G','2SA':'N','2SI':'X','2ST':'T',
    '2TL':'T','2TY':'Y','2VA':'V','2XA':'C','32S':'X',
    '32T':'X','3AH':'H','3AR':'X','3CF':'F','3DA':'A',
    '3DR':'N','3GA':'A','3MD':'D','3ME':'U','3NF':'Y',
    '3QN':'K','3TY':'X','3XH':'G','4AC':'N','4BF':'Y',
    '4CF':'F','4CY':'M','4DP':'W','4F3':'GYG','4FB':'P',
    '4FW':'W','4HT':'W','4IN':'W','4MF':'N','4MM':'X',
    '4OC':'C','4PC':'C','4PD':'C','4PE':'C','4PH':'F',
    '4SC':'C','4SU':'U','4TA':'N','4U7':'A','56A':'H',
    '5AA':'A','5AB':'A','5AT':'T','5BU':'U','5CG':'G',
    '5CM':'C','5CS':'C','5FA':'A','5FC':'C','5FU':'U',
    '5HP':'E','5HT':'T','5HU':'U','5IC':'C','5IT':'T',
    '5IU':'U','5MC':'C','5MD':'N','5MU':'U','5NC':'C',
    '5PC':'C','5PY':'T','5SE':'U','5ZA':'TWG','64T':'T',
    '6CL':'K','6CT':'T','6CW':'W','6HA':'A','6HC':'C',
    '6HG':'G','6HN':'K','6HT':'T','6IA':'A','6MA':'A',
    '6MC':'A','6MI':'N','6MT':'A','6MZ':'N','6OG':'G',
    '70U':'U','7DA':'A','7GU':'G','7JA':'I','7MG':'G',
    '8AN':'A','8FG':'G','8MG':'G','8OG':'G','9NE':'E',
    '9NF':'F','9NR':'R','9NV':'V','A  ':'A','A1P':'N',
    'A23':'A','A2L':'A','A2M':'A','A34':'A','A35':'A',
    'A38':'A','A39':'A','A3A':'A','A3P':'A','A40':'A',
    'A43':'A','A44':'A','A47':'A','A5L':'A','A5M':'C',
    'A5N':'N','A5O':'A','A66':'X','AA3':'A','AA4':'A',
    'AAR':'R','AB7':'X','ABA':'A','ABR':'A','ABS':'A',
    'ABT':'N','ACB':'D','ACL':'R','AD2':'A','ADD':'X',
    'ADX':'N','AEA':'X','AEI':'D','AET':'A','AFA':'N',
    'AFF':'N','AFG':'G','AGM':'R','AGT':'C','AHB':'N',
    'AHH':'X','AHO':'A','AHP':'A','AHS':'X','AHT':'X',
    'AIB':'A','AKL':'D','AKZ':'D','ALA':'A','ALC':'A',
    'ALM':'A','ALN':'A','ALO':'T','ALQ':'X','ALS':'A',
    'ALT':'A','ALV':'A','ALY':'K','AN8':'A','AP7':'A',
    'APE':'X','APH':'A','API':'K','APK':'K','APM':'X',
    'APP':'X','AR2':'R','AR4':'E','AR7':'R','ARG':'R',
    'ARM':'R','ARO':'R','ARV':'X','AS ':'A','AS2':'D',
    'AS9':'X','ASA':'D','ASB':'D','ASI':'D','ASK':'D',
    'ASL':'D','ASM':'X','ASN':'N','ASP':'D','ASQ':'D',
    'ASU':'N','ASX':'B','ATD':'T','ATL':'T','ATM':'T',
    'AVC':'A','AVN':'X','AYA':'A','AYG':'AYG','AZK':'K',
    'AZS':'S','AZY':'Y','B1F':'F','B1P':'N','B2A':'A',
    'B2F':'F','B2I':'I','B2V':'V','B3A':'A','B3D':'D',
    'B3E':'E','B3K':'K','B3L':'X','B3M':'X','B3Q':'X',
    'B3S':'S','B3T':'X','B3U':'H','B3X':'N','B3Y':'Y',
    'BB6':'C','BB7':'C','BB8':'F','BB9':'C','BBC':'C',
    'BCS':'C','BE2':'X','BFD':'D','BG1':'S','BGM':'G',
    'BH2':'D','BHD':'D','BIF':'F','BIL':'X','BIU':'I',
    'BJH':'X','BLE':'L','BLY':'K','BMP':'N','BMT':'T',
    'BNN':'F','BNO':'X','BOE':'T','BOR':'R','BPE':'C',
    'BRU':'U','BSE':'S','BT5':'N','BTA':'L','BTC':'C',
    'BTR':'W','BUC':'C','BUG':'V','BVP':'U','BZG':'N',
    'C  ':'C','C12':'TYG','C1X':'K','C25':'C','C2L':'C',
    'C2S':'C','C31':'C','C32':'C','C34':'C','C36':'C',
    'C37':'C','C38':'C','C3Y':'C','C42':'C','C43':'C',
    'C45':'C','C46':'C','C49':'C','C4R':'C','C4S':'C',
    'C5C':'C','C66':'X','C6C':'C','C99':'TFG','CAF':'C',
    'CAL':'X','CAR':'C','CAS':'C','CAV':'X','CAY':'C',
    'CB2':'C','CBR':'C','CBV':'C','CCC':'C','CCL':'K',
    'CCS':'C','CCY':'CYG','CDE':'X','CDV':'X','CDW':'C',
    'CEA':'C','CFL':'C','CFY':'FCYG','CG1':'G','CGA':'E',
    'CGU':'E','CH ':'C','CH6':'MYG','CH7':'KYG','CHF':'X',
    'CHG':'X','CHP':'G','CHS':'X','CIR':'R','CJO':'GYG',
    'CLE':'L','CLG':'K','CLH':'K','CLV':'AFG','CM0':'N',
    'CME':'C','CMH':'C','CML':'C','CMR':'C','CMT':'C',
    'CNU':'U','CP1':'C','CPC':'X','CPI':'X','CQR':'GYG',
    'CR0':'TLG','CR2':'GYG','CR5':'G','CR7':'KYG','CR8':'HYG',
    'CRF':'TWG','CRG':'THG','CRK':'MYG','CRO':'GYG','CRQ':'QYG',
    'CRU':'EYG','CRW':'ASG','CRX':'ASG','CS0':'C','CS1':'C',
    'CS3':'C','CS4':'C','CS8':'N','CSA':'C','CSB':'C',
    'CSD':'C','CSE':'C','CSF':'C','CSH':'SHG','CSI':'G',
    'CSJ':'C','CSL':'C','CSO':'C','CSP':'C','CSR':'C',
    'CSS':'C','CSU':'C','CSW':'C','CSX':'C','CSY':'SYG',
    'CSZ':'C','CTE':'W','CTG':'T','CTH':'T','CUC':'X',
    'CWR':'S','CXM':'M','CY0':'C','CY1':'C','CY3':'C',
    'CY4':'C','CYA':'C','CYD':'C','CYF':'C','CYG':'C',
    'CYJ':'X','CYM':'C','CYQ':'C','CYR':'C','CYS':'C',
    'CZ2':'C','CZO':'GYG','CZZ':'C','D11':'T','D1P':'N',
    'D3 ':'N','D33':'N','D3P':'G','D3T':'T','D4M':'T',
    'D4P':'X','DA ':'A','DA2':'X','DAB':'A','DAH':'F',
    'DAL':'A','DAR':'R','DAS':'D','DBB':'T','DBM':'N',
    'DBS':'S','DBU':'T','DBY':'Y','DBZ':'A','DC ':'C',
    'DC2':'C','DCG':'G','DCI':'X','DCL':'X','DCT':'C',
    'DCY':'C','DDE':'H','DDG':'G','DDN':'U','DDX':'N',
    'DFC':'C','DFG':'G','DFI':'X','DFO':'X','DFT':'N',
    'DG ':'G','DGH':'G','DGI':'G','DGL':'E','DGN':'Q',
    'DHA':'S','DHI':'H','DHL':'X','DHN':'V','DHP':'X',
    'DHU':'U','DHV':'V','DI ':'I','DIL':'I','DIR':'R',
    'DIV':'V','DLE':'L','DLS':'K','DLY':'K','DM0':'K',
    'DMH':'N','DMK':'D','DMT':'X','DN ':'N','DNE':'L',
    'DNG':'L','DNL':'K','DNM':'L','DNP':'A','DNR':'C',
    'DNS':'K','DOA':'X','DOC':'C','DOH':'D','DON':'L',
    'DPB':'T','DPH':'F','DPL':'P','DPP':'A','DPQ':'Y',
    'DPR':'P','DPY':'N','DRM':'U','DRP':'N','DRT':'T',
    'DRZ':'N','DSE':'S','DSG':'N','DSN':'S','DSP':'D',
    'DT ':'T','DTH':'T','DTR':'W','DTY':'Y','DU ':'U',
    'DVA':'V','DXD':'N','DXN':'N','DYG':'DYG','DYS':'C',
    'DZM':'A','E  ':'A','E1X':'A','ECC':'Q','EDA':'A',
    'EFC':'C','EHP':'F','EIT':'T','ENP':'N','ESB':'Y',
    'ESC':'M','EXB':'X','EXY':'L','EY5':'N','EYS':'X',
    'F2F':'F','FA2':'A','FA5':'N','FAG':'N','FAI':'N',
    'FB5':'A','FB6':'A','FCL':'F','FFD':'N','FGA':'E',
    'FGL':'G','FGP':'S','FHL':'X','FHO':'K','FHU':'U',
    'FLA':'A','FLE':'L','FLT':'Y','FME':'M','FMG':'G',
    'FMU':'N','FOE':'C','FOX':'G','FP9':'P','FPA':'F',
    'FRD':'X','FT6':'W','FTR':'W','FTY':'Y','FVA':'V',
    'FZN':'K','G  ':'G','G25':'G','G2L':'G','G2S':'G',
    'G31':'G','G32':'G','G33':'G','G36':'G','G38':'G',
    'G42':'G','G46':'G','G47':'G','G48':'G','G49':'G',
    'G4P':'N','G7M':'G','GAO':'G','GAU':'E','GCK':'C',
    'GCM':'X','GDP':'G','GDR':'G','GFL':'G','GGL':'E',
    'GH3':'G','GHG':'Q','GHP':'G','GL3':'G','GLH':'Q',
    'GLJ':'E','GLK':'E','GLM':'X','GLN':'Q','GLQ':'E',
    'GLU':'E','GLX':'Z','GLY':'G','GLZ':'G','GMA':'E',
    'GMS':'G','GMU':'U','GN7':'G','GND':'X','GNE':'N',
    'GOM':'G','GPL':'K','GS ':'G','GSC':'G','GSR':'G',
    'GSS':'G','GSU':'E','GT9':'C','GTP':'G','GVL':'X',
    'GYC':'CYG','GYS':'SYG','H2U':'U','H5M':'P','HAC':'A',
    'HAR':'R','HBN':'H','HCS':'X','HDP':'U','HEU':'U',
    'HFA':'X','HGL':'X','HHI':'H','HHK':'AK','HIA':'H',
    'HIC':'H','HIP':'H','HIQ':'H','HIS':'H','HL2':'L',
    'HLU':'L','HMR':'R','HOL':'N','HPC':'F','HPE':'F',
    'HPH':'F','HPQ':'F','HQA':'A','HRG':'R','HRP':'W',
    'HS8':'H','HS9':'H','HSE':'S','HSL':'S','HSO':'H',
    'HTI':'C','HTN':'N','HTR':'W','HV5':'A','HVA':'V',
    'HY3':'P','HYP':'P','HZP':'P','I  ':'I','I2M':'I',
    'I58':'K','I5C':'C','IAM':'A','IAR':'R','IAS':'D',
    'IC ':'C','IEL':'K','IEY':'HYG','IG ':'G','IGL':'G',
    'IGU':'G','IIC':'SHG','IIL':'I','ILE':'I','ILG':'E',
    'ILX':'I','IMC':'C','IML':'I','IOY':'F','IPG':'G',
    'IPN':'N','IRN':'N','IT1':'K','IU ':'U','IYR':'Y',
    'IYT':'T','IZO':'M','JJJ':'C','JJK':'C','JJL':'C',
    'JW5':'N','K1R':'C','KAG':'G','KCX':'K','KGC':'K',
    'KNB':'A','KOR':'M','KPI':'K','KST':'K','KYQ':'K',
    'L2A':'X','LA2':'K','LAA':'D','LAL':'A','LBY':'K',
    'LC ':'C','LCA':'A','LCC':'N','LCG':'G','LCH':'N',
    'LCK':'K','LCX':'K','LDH':'K','LED':'L','LEF':'L',
    'LEH':'L','LEI':'V','LEM':'L','LEN':'L','LET':'X',
    'LEU':'L','LEX':'L','LG ':'G','LGP':'G','LHC':'X',
    'LHU':'U','LKC':'N','LLP':'K','LLY':'K','LME':'E',
    'LMF':'K','LMQ':'Q','LMS':'N','LP6':'K','LPD':'P',
    'LPG':'G','LPL':'X','LPS':'S','LSO':'X','LTA':'X',
    'LTR':'W','LVG':'G','LVN':'V','LYF':'K','LYK':'K',
    'LYM':'K','LYN':'K','LYR':'K','LYS':'K','LYX':'K',
    'LYZ':'K','M0H':'C','M1G':'G','M2G':'G','M2L':'K',
    'M2S':'M','M30':'G','M3L':'K','M5M':'C','MA ':'A',
    'MA6':'A','MA7':'A','MAA':'A','MAD':'A','MAI':'R',
    'MBQ':'Y','MBZ':'N','MC1':'S','MCG':'X','MCL':'K',
    'MCS':'C','MCY':'C','MD3':'C','MD6':'G','MDH':'X',
    'MDO':'ASG','MDR':'N','MEA':'F','MED':'M','MEG':'E',
    'MEN':'N','MEP':'U','MEQ':'Q','MET':'M','MEU':'G',
    'MF3':'X','MFC':'GYG','MG1':'G','MGG':'R','MGN':'Q',
    'MGQ':'A','MGV':'G','MGY':'G','MHL':'L','MHO':'M',
    'MHS':'H','MIA':'A','MIS':'S','MK8':'L','ML3':'K',
    'MLE':'L','MLL':'L','MLY':'K','MLZ':'K','MME':'M',
    'MMO':'R','MMT':'T','MND':'N','MNL':'L','MNU':'U',
    'MNV':'V','MOD':'X','MP8':'P','MPH':'X','MPJ':'X',
    'MPQ':'G','MRG':'G','MSA':'G','MSE':'M','MSL':'M',
    'MSO':'M','MSP':'X','MT2':'M','MTR':'T','MTU':'A',
    'MTY':'Y','MVA':'V','N  ':'N','N10':'S','N2C':'X',
    'N5I':'N','N5M':'C','N6G':'G','N7P':'P','NA8':'A',
    'NAL':'A','NAM':'A','NB8':'N','NBQ':'Y','NC1':'S',
    'NCB':'A','NCX':'N','NCY':'X','NDF':'F','NDN':'U',
    'NEM':'H','NEP':'H','NF2':'N','NFA':'F','NHL':'E',
    'NIT':'X','NIY':'Y','NLE':'L','NLN':'L','NLO':'L',
    'NLP':'L','NLQ':'Q','NMC':'G','NMM':'R','NMS':'T',
    'NMT':'T','NNH':'R','NP3':'N','NPH':'C','NPI':'A',
    'NRP':'LYG','NRQ':'MYG','NSK':'X','NTY':'Y','NVA':'V',
    'NYC':'TWG','NYG':'NYG','NYM':'N','NYS':'C','NZH':'H',
    'O12':'X','O2C':'N','O2G':'G','OAD':'N','OAS':'S',
    'OBF':'X','OBS':'X','OCS':'C','OCY':'C','ODP':'N',
    'OHI':'H','OHS':'D','OIC':'X','OIP':'I','OLE':'X',
    'OLT':'T','OLZ':'S','OMC':'C','OMG':'G','OMT':'M',
    'OMU':'U','ONE':'U','ONH':'A','ONL':'X','OPR':'R',
    'ORN':'A','ORQ':'R','OSE':'S','OTB':'X','OTH':'T',
    'OTY':'Y','OXX':'D','P  ':'G','P1L':'C','P1P':'N',
    'P2T':'T','P2U':'U','P2Y':'P','P5P':'A','PAQ':'Y',
    'PAS':'D','PAT':'W','PAU':'A','PBB':'C','PBF':'F',
    'PBT':'N','PCA':'E','PCC':'P','PCE':'X','PCS':'F',
    'PDL':'X','PDU':'U','PEC':'C','PF5':'F','PFF':'F',
    'PFX':'X','PG1':'S','PG7':'G','PG9':'G','PGL':'X',
    'PGN':'G','PGP':'G','PGY':'G','PHA':'F','PHD':'D',
    'PHE':'F','PHI':'F','PHL':'F','PHM':'F','PIA':'AYG',
    'PIV':'X','PLE':'L','PM3':'F','PMT':'C','POM':'P',
    'PPN':'F','PPU':'A','PPW':'G','PQ1':'N','PR3':'C',
    'PR5':'A','PR9':'P','PRN':'A','PRO':'P','PRS':'P',
    'PSA':'F','PSH':'H','PST':'T','PSU':'U','PSW':'C',
    'PTA':'X','PTH':'Y','PTM':'Y','PTR':'Y','PU ':'A',
    'PUY':'N','PVH':'H','PVL':'X','PYA':'A','PYO':'U',
    'PYX':'C','PYY':'N','QLG':'QLG','QMM':'Q','QPA':'C',
    'QPH':'F','QUO':'G','R  ':'A','R1A':'C','R4K':'W',
    'RC7':'HYG','RE0':'W','RE3':'W','RIA':'A','RMP':'A',
    'RON':'X','RT ':'T','RTP':'N','S1H':'S','S2C':'C',
    'S2D':'A','S2M':'T','S2P':'A','S4A':'A','S4C':'C',
    'S4G':'G','S4U':'U','S6G':'G','SAC':'S','SAH':'C',
    'SAR':'G','SBL':'S','SC ':'C','SCH':'C','SCS':'C',
    'SCY':'C','SD2':'X','SDG':'G','SDP':'S','SEB':'S',
    'SEC':'A','SEG':'A','SEL':'S','SEM':'S','SEN':'S',
    'SEP':'S','SER':'S','SET':'S','SGB':'S','SHC':'C',
    'SHP':'G','SHR':'K','SIB':'C','SIC':'DC','SLA':'P',
    'SLR':'P','SLZ':'K','SMC':'C','SME':'M','SMF':'F',
    'SMP':'A','SMT':'T','SNC':'C','SNN':'N','SOC':'C',
    'SOS':'N','SOY':'S','SPT':'T','SRA':'A','SSU':'U',
    'STY':'Y','SUB':'X','SUI':'DG','SUN':'S','SUR':'U',
    'SVA':'S','SVV':'S','SVW':'S','SVX':'S','SVY':'S',
    'SVZ':'X','SWG':'SWG','SYS':'C','T  ':'T','T11':'F',
    'T23':'T','T2S':'T','T2T':'N','T31':'U','T32':'T',
    'T36':'T','T37':'T','T38':'T','T39':'T','T3P':'T',
    'T41':'T','T48':'T','T49':'T','T4S':'T','T5O':'U',
    'T5S':'T','T66':'X','T6A':'A','TA3':'T','TA4':'X',
    'TAF':'T','TAL':'N','TAV':'D','TBG':'V','TBM':'T',
    'TC1':'C','TCP':'T','TCQ':'Y','TCR':'W','TCY':'A',
    'TDD':'L','TDY':'T','TFE':'T','TFO':'A','TFQ':'F',
    'TFT':'T','TGP':'G','TH6':'T','THC':'T','THO':'X',
    'THR':'T','THX':'N','THZ':'R','TIH':'A','TLB':'N',
    'TLC':'T','TLN':'U','TMB':'T','TMD':'T','TNB':'C',
    'TNR':'S','TOX':'W','TP1':'T','TPC':'C','TPG':'G',
    'TPH':'X','TPL':'W','TPO':'T','TPQ':'Y','TQI':'W',
    'TQQ':'W','TRF':'W','TRG':'K','TRN':'W','TRO':'W',
    'TRP':'W','TRQ':'W','TRW':'W','TRX':'W','TS ':'N',
    'TST':'X','TT ':'N','TTD':'T','TTI':'U','TTM':'T',
    'TTQ':'W','TTS':'Y','TY1':'Y','TY2':'Y','TY3':'Y',
    'TY5':'Y','TYB':'Y','TYI':'Y','TYJ':'Y','TYN':'Y',
    'TYO':'Y','TYQ':'Y','TYR':'Y','TYS':'Y','TYT':'Y',
    'TYU':'N','TYW':'Y','TYX':'X','TYY':'Y','TZB':'X',
    'TZO':'X','U  ':'U','U25':'U','U2L':'U','U2N':'U',
    'U2P':'U','U31':'U','U33':'U','U34':'U','U36':'U',
    'U37':'U','U8U':'U','UAR':'U','UCL':'U','UD5':'U',
    'UDP':'N','UFP':'N','UFR':'U','UFT':'U','UMA':'A',
    'UMP':'U','UMS':'U','UN1':'X','UN2':'X','UNK':'X',
    'UR3':'U','URD':'U','US1':'U','US2':'U','US3':'T',
    'US5':'U','USM':'U','VAD':'V','VAF':'V','VAL':'V',
    'VB1':'K','VDL':'X','VLL':'X','VLM':'X','VMS':'X',
    'VOL':'X','WCR':'GYG','X  ':'G','X2W':'E','X4A':'N',
    'X9Q':'AFG','XAD':'A','XAE':'N','XAL':'A','XAR':'N',
    'XCL':'C','XCN':'C','XCP':'X','XCR':'C','XCS':'N',
    'XCT':'C','XCY':'C','XGA':'N','XGL':'G','XGR':'G',
    'XGU':'G','XPR':'P','XSN':'N','XTH':'T','XTL':'T',
    'XTR':'T','XTS':'G','XTY':'N','XUA':'A','XUG':'G',
    'XX1':'K','XXY':'THG','XYG':'DYG','Y  ':'A','YCM':'C',
    'YG ':'G','YOF':'Y','YRR':'N','YYG':'G','Z  ':'C',
    'Z01':'A','ZAD':'A','ZAL':'A','ZBC':'C','ZBU':'U',
    'ZCL':'F','ZCY':'C','ZDU':'U','ZFB':'X','ZGU':'G',
    'ZHP':'N','ZTH':'T','ZU0':'T','ZZJ':'A', 'DT':'T', 'DA':'A', 'DC':'C', 'DU': 'U', 'DG':'G' }

