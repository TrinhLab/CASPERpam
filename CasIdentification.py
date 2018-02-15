"""For taking the information from Makarova 2015 and gaining the cas identities of all species listed."""

f = open("/Users/brianmendoza/Desktop/CasIDs.txt")
species_and_type = dict()
for line in f:
    # only those that start with a '#' will have the genome name and the type/subtype included
    if line[0] == '#':
        mylist = line.split('\t')
        species_list = mylist[3].split("_")
        species = species_list[0] + "_" + species_list[1]
        # To check to see if the species has already been seen
        if species in species_and_type:
            species_and_type[species] += ";" + mylist[1]
        else:
            species_and_type[species] = mylist[1]

        # determine which types the species has
        cas_type = species_and_type[species]
        is_forward = False
        is_back = False
        if cas_type.find('CAS-I-') != -1 or cas_type.find('CAS-III-') != -1 or cas_type.find('CAS-IV-') != -1:
            is_forward = True
        elif cas_type.find('CAS-II-') or cas_type.find('Cas-V-'):
            is_back = True
        else:
            is_back = True
            is_forward = True
        # assign pam_location
        if is_back and is_forward:
            pam_location = 'both'
        elif is_back:
            pam_location = 'back'
        else:
            pam_location = 'forward'
f.close()
for item in species_and_type:
    print(item + ";" + species_and_type[item])
