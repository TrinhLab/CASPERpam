"""This is used to get the confidence score for every species that have a blast hit (2,990)."""

import ProcessPAMbySpecies

# A dictionary lists, stores the Cas-type information at the genus level
genus_dict = dict()

# a dictionary by species that keeps track of the score for each species
score_tracker = dict()


# uses the information from Makarova et al. 2015. This also generates the Cm/Ct score
def import_makarova_data():
    f = open('/Users/brianmendoza/Dropbox/Code_Transfer/cas_type_array_Makarova.txt')
    # Transfer the information into genus_dict
    for line in f:
        species_list = line.split("\t")
        genus = species_list[17]
        species = species_list[18][:-1]
        if genus in genus_dict:
            genus_dict[genus][species] = list()
        else:
            genus_dict[genus] = dict()
            genus_dict[genus][species] = list()
        for i in range(16):
            genus_dict[genus][species].append(int(species_list[i]))
    f.close()


# gets the information from pamdiscovery_complete.
def import_species_hit_data():
    ProcessPAMbySpecies.build_hit_database(scoring=True)
    final_dict = ProcessPAMbySpecies.final_dict
    for species in final_dict:
        hit_genomes = set()
        spacers = set()
        total_hits = 0
        num_assemblies  = len(final_dict[species])
        for assembly in final_dict[species]:
            for hit in final_dict[species][assembly]:
                hit_genomes.add(hit[1])
                spacers.add(hit[0])
                total_hits += 1
        # initialize the score tracker
        score_tracker[species] = [total_hits, num_assemblies, len(hit_genomes), len(spacers)]


def confidence_scores():
    for species in score_tracker:

        # This section of code gets the Cm/CT score for the similarity of Cas Type
        if species[0] in genus_dict:
            if species[1] in genus_dict[species[0]]:
                caslist = genus_dict[species[0]][species[1]]
                m = max(caslist)
                T = sum(caslist)
                print(species, m, T)
                cmct = m/T
            else:
                totcaslist = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
                for spec in genus_dict[species[0]]:
                    for i in range(len(totcaslist)):
                        totcaslist[i] += genus_dict[species[0]][spec][i]
                m = max(totcaslist)
                T = sum(totcaslist)
                cmct = m/(T*2)
        else:
            cmct = 0.125
        score_tracker[species].append(cmct)

        # This section of code consolidates all the scores
        # going to leave this to excel for now so we can see the granularity of each score


def output():
    f = open('confidence_scores.txt','w')
    for species in score_tracker:
        result_string = str()
        for item in score_tracker[species]:
            result_string += "," + str(item)
        f.write(species[0] + "_" + species[1] + result_string + '\n')
    f.close()


# Execution:
import_species_hit_data()
import_makarova_data()
confidence_scores()
output()
