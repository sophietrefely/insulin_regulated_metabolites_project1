import requests
response = requests.get('http://reactomews.oicr.on.ca:8080/ReactomeRESTfulAPI/RESTfulWS/pathwayParticipants/5934172')#example for PFKFB3
response.json()
for result in response.json():
    if result[u'schemaClass'] == 'SimpleEntity':
        print(result['displayName'])


reactome_id_list = [5934172, 3099]

import requests
metabolite_dict = {}
for item in reactome_id_list:
    reactome_id = str(item)
    url = 'http://reactomews.oicr.on.ca:8080/ReactomeRESTfulAPI/RESTfulWS/pathwayParticipants/' + reactome_id
    response = requests.get(url)
    response.json() #convert to JSON
    for result in response.json():
        if result[u'schemaClass'] == 'SimpleEntity': #response for all compounds is 'SimpleEntity', whereas proteins with interactions are 'DefinedSet'.
            print(result['displayName'])
            print(result['dbId'])
            metabolite = result['displayName'] #store metabolite and its reactome ID
            metabolite_dict[metabolite] = metabolite_dict.get(metabolite, 0) + 1  #use the dictionary to count amount of occurances for each metabolite.
print(metabolite_dict)
