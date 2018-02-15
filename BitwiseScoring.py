"""This file generates the bitwise information for each species and obtains PAM information for every species"""
import math, numpy


all_pams = dict()
species_type_dict = dict()

#Output:
final_pams_scores = dict()


# IMPORT FUNCTIONS
def import_pams():
    f = open("hitsbyspecies4.txt")
    for line in f:
        itemlist = line[:-2].split(",")
        species = itemlist[1]
        flanks = list()
        for i in range(2,len(itemlist)):
            flanks.append(itemlist[i])
        all_pams[species] = flanks
    f.close()


def import_cas_categories():
    f = open("/Users/brianmendoza/Desktop/total_cas_categories.txt")
    for line in f:
       myline = line.split("\t")
       species_type_dict[myline[0]] = myline[1]
    f.close()

# ----------------------------END IMPORT FUNCTIONS --------------------------- #


# SCORING FUNCTIONS FOR BITWISE INFORMATION
def pams_and_confidence():
    for item in species_type_dict:
        # figure out whether to use back or front flanks
        if species_type_dict[item].find("-II-") == -1:
            flank_type = "_front"
            five_prime = True
        else:
            flank_type = "_back"
            five_prime = False
        cur_species = item + flank_type
        metrics = get_metrics(all_pams[cur_species], five_prime)
        final_pams_scores[cur_species] = metrics


def get_metrics(flank_list, five_prime):
    # species-wide variables
    N = len(flank_list)
    ni = len(flank_list[0])
    en = 1.039721/N
    # i specific variables
    fai_matrix = list()
    Ri_list = list()
    for i in range(0,ni):
        freq = {"A":0,"T":0,"C":0,"G":0}
        for flank in flank_list:
            freq[flank[i]] += 1
        # get relative frequencies for each nucleotide
        for item in freq:
            freq[item] = freq[item]/N
        fai_matrix.append(freq)
        Ri_list.append(get_ri(freq,en))
    bby_pam = get_pam(Ri_list,fai_matrix,five_prime)
    return bby_pam


def get_ri(freq_table,en):
    hi_score = 0
    for item in freq_table:
        # to make sure the log doesn't see a zero:
        if freq_table[item] != 0:
            hi_score += freq_table[item]*math.log(freq_table[item],2)
    ri_score = 2+hi_score-en
    return ri_score


def get_pam(Ris,fais,is_five_prime):
    pam_nt_and_pos = dict()
    risdev = list()
    # Find the significant positions and put them into the sig_pos container
    sig_pos = list()
    Ri_mean = numpy.mean(Ris)
    Ri_stdev = numpy.std(Ris)
    for i in range(len(Ris)):
        if Ris[i] > Ri_mean+Ri_stdev/2:
            # Check if it is close to sequence PAM is unlikely past first 5 nucleotides
            if is_five_prime:
                if i < len(Ris)-5:
                    if Ris[i] > Ri_mean + Ri_stdev:
                        sig_pos.append(i)
                        risdev.append(i)
                else:
                    sig_pos.append(i)
                    risdev.append(Ris[i] - Ri_mean)
            else:
                if i > 5:
                    if Ris[i] > Ri_mean + Ri_stdev:
                        sig_pos.append(i)
                        risdev.append(i)
                else:
                    sig_pos.append(i)
                    risdev.append(Ris[i]-Ri_mean)
    # check to see if there were any significant sequences at all:
    if not sig_pos:
        return "No PAM identified."
    for pos in sig_pos:
        for nt in fais[pos]:
            # single consensus nucleotide
            if fais[pos][nt] > 0.5:
                pam_nt_and_pos[pos] = nt
            # possible degenerate code
            elif fais[pos][nt] > 0.25:
                if pos in pam_nt_and_pos:
                    pam_nt_and_pos[pos] += nt
                else:
                    pam_nt_and_pos[pos] = nt

    # find the start and end of the pam depending on the type:
    if is_five_prime:
        pam_start = min(pam_nt_and_pos.keys())
        pam_end = len(Ris)
    else:
        pam_start = 0
        pam_end = max(pam_nt_and_pos.keys())+1

    # Make the PAM sequence by putting N's in the places that are not significant
    PAM = str()
    for i in range(pam_start,pam_end):
        if i in pam_nt_and_pos:
            if len(pam_nt_and_pos[i]) > 1:
                PAM += degenerate_nucleotides(pam_nt_and_pos[i])
            else:
                PAM += pam_nt_and_pos[i]
        else:
            PAM += "N"
    return PAM, risdev


# --------------------------END IMPORT FUNCTIONS------------------------------------------- #


def export_by_species():
    f = open("PAMscores.txt", "w")
    for species in final_pams_scores:
        bby_sum = 0.0
        totris = len(final_pams_scores[species][1])
        for item in final_pams_scores[species][1]:
            if type(item) is not str:
                bby_sum += item
        f.write(species + "," + final_pams_scores[species][0] + "," + str(bby_sum/totris) + "\n")
    f.close()

# ----------------------- HELPER FUNCTIONS ------------- #


def degenerate_nucleotides(nts):
    if nts == "AG":
        return "R"
    elif nts == "TC":
        return "Y"
    elif nts == "CG":
        return "S"
    elif nts == "AT":
        return "W"
    elif nts == "TG":
        return "K"
    elif nts == "AC":
        return "M"
    elif nts.find("A") == -1:
        return "B"
    elif nts.find("T") == -1:
        return "V"
    elif nts.find("C") == -1:
        return "D"
    elif nts.find("G") == -1:
        return "H"
    else:
        return "N"


# ------- EXECUTION ----- #
import_pams()
import_cas_categories()
pams_and_confidence()
export_by_species()
