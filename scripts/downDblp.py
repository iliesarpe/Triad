def download(url, filename):
	import os
	import requests
	from tqdm import tqdm

	print("Downloading", filename)
	if not os.path.isfile(filename):
		print("file is not there:", filename)
		exit(0)
		r = requests.get(url)
		total_size_in_bytes= int(r.headers.get('content-length', 0))
		print("downloading", total_size_in_bytes, "bytes")
		progress = tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True)
		with open(filename, 'wb') as fd:
			for chunk in r.iter_content(chunk_size=4096):
				progress.update(len(chunk))
				fd.write(chunk)
		progress.close()

def dataset_dblp():
	import itertools
	import xml.dom.pulldom as pdom
	import gzip
	import json
	from tqdm import tqdm

	url = 'https://dblp.org/xml/release/dblp-2024-08-01.xml.gz'
	local_xml = 'data/dblp-2024-08-01.xml.gz'
	local_raw_csv = 'data/dblp-2024-08-01.raw.csv'
	local_mapping_json = 'data/dblp-2024-08-01.raw-authorsMap.json'
	download(url, local_xml)

	names = dict()
	def remap_name(name):
		if name in names:
			return names[name]
		else:
			identifier = len(names)
			names[name] = identifier
			return identifier

	output = open(local_raw_csv, 'w')
	active_items = ['article', 'inproceedings', 'incollection', 'book']
	i = 100000
	authors_venues = {}
	existing_edges = {}
	with tqdm.wrapattr(gzip.open(local_xml), 'read') as fp:
		doc = pdom.parse(fp)
		for event, node in doc:
			if event == pdom.START_ELEMENT and node.tagName in active_items:
				doc.expandNode(node)
				year = node.getElementsByTagName('year')
				venue = None
				if event.find('booktitle'):
					if len(node.getElementsByTagName('booktitle')) > 0:
						for ch in node.getElementsByTagName('booktitle').item(0).childNodes:
							if venue is None:
								venue = ch.data
							else:
								venue += ch.data + " "
						#	print(ch.data, end='')
						#venue = node.getElementsByTagName('booktitle').item(0).firstChild.data
				elif event.find('journal'):
					if len(node.getElementsByTagName('journal')) > 0:
						for ch in node.getElementsByTagName('booktitle').item(0).childNodes:
							if venue is None:
								venue = ch.data
							else:
								venue += ch.data + " "
						#venue = node.getElementsByTagName('journal').item(0).data
				#i -=1
				#if i ==0:
				#	break
				if venue is not None:
					venue = venue.replace(",", "")
					authors = []
					for auth in node.getElementsByTagName('author'):
						if auth.firstChild.data not in authors_venues.keys():
							authors_venues[str(auth.firstChild.data)] = {venue: 1}
						else:
							if venue not in authors_venues[str(auth.firstChild.data)].keys():
								authors_venues[str(auth.firstChild.data)][venue] = 1
							else:
								authors_venues[str(auth.firstChild.data)][venue] += 1
					if len(node.getElementsByTagName('author')) > 1:
						authors = [
							remap_name(a.firstChild.data)
							for a in node.getElementsByTagName('author')
						]
					if len(authors) > 1:
						authors.sort()
					for a, b in itertools.combinations(authors, 2):
						src = min(a,b)
						dst = max(a,b)
						if src not in existing_edges.keys():
							existing_edges[src] = set([dst])
							print(f"{src} {dst}", file=output)
						elif dst not in existing_edges[src]:
							existing_edges[src].add(dst)
							print(f"{src} {dst}", file=output)
	output.close()
	dict_mapping = {"authormapping": names, "authorsVenues": authors_venues}
	json_obj = json.dumps(dict_mapping)
	with open(local_mapping_json, "w") as outj:
		outj.write(json_obj)

def dataset_temporal_dblp():
	import itertools
	import xml.dom.pulldom as pdom
	import gzip
	import json
	from tqdm import tqdm

	url = 'https://dblp.org/xml/release/dblp-2024-08-01.xml.gz'
	local_xml = 'data/dblp-2024-08-01.xml.gz'
	local_raw_csv = 'data/dblp-2024-08-01.raw.csv'
	local_mapping_json = 'data/dblp-2024-08-01-SnapshotsYears.json'

	years = [[a,a+5] for a in range(1970, 2024, 5)]
	year_map = {}
	for i,year_range in enumerate(years):
		for y in range(year_range[0], year_range[1]):
			year_map[y] = i
	#for year in year_map.keys():
	#	print(year, year_map[year])

	authors_venues_years = dict({})
	existing_edges_years = dict({})
	existing_edges_years_alg = dict({})
	for i in range(len(years)):
		authors_venues_years[i] = dict({})
		existing_edges_years[i] = dict({})
		existing_edges_years_alg[i] = dict({})
		
	#exit(0)

	#download(url, local_xml)

	names = dict()
	def remap_name(name):
		if name in names:
			return names[name]
		else:
			identifier = len(names)
			names[name] = identifier
			return identifier

	#output = open(local_raw_csv, 'w')
	active_items = ['article', 'inproceedings', 'incollection', 'book']
	i = 100000
	#authors_venues = {}

	with tqdm.wrapattr(gzip.open(local_xml), 'read') as fp:
		doc = pdom.parse(fp)
		for event, node in doc:
			if event == pdom.START_ELEMENT and node.tagName in active_items:
				doc.expandNode(node)
				year = node.getElementsByTagName('year')
				if year is not None and len(year) > 0:
					year = int(year[0].firstChild.data)
				else:
					continue
				if year not in year_map.keys():
					continue
				year_key = year_map[year]
				venue = None
				if event.find('booktitle'):
					if len(node.getElementsByTagName('booktitle')) > 0:
						for ch in node.getElementsByTagName('booktitle').item(0).childNodes:
							if venue is None:
								venue = ch.data
							else:
								venue += ch.data + " "
						#	print(ch.data, end='')
						#venue = node.getElementsByTagName('booktitle').item(0).firstChild.data
				elif event.find('journal'):
					if len(node.getElementsByTagName('journal')) > 0:
						for ch in node.getElementsByTagName('booktitle').item(0).childNodes:
							if venue is None:
								venue = ch.data
							else:
								venue += ch.data + " "
						#venue = node.getElementsByTagName('journal').item(0).data
				#i -=1
				#if i ==0:
				#	break
				if venue is not None:
					venue = venue.replace(",", "")
					authors = []
					for auth in node.getElementsByTagName('author'):
						if auth.firstChild.data not in authors_venues_years[year_key].keys():
							authors_venues_years[year_key][str(auth.firstChild.data)] = {venue: 1}
						else:
							if venue not in authors_venues_years[year_key][str(auth.firstChild.data)].keys():
								authors_venues_years[year_key][str(auth.firstChild.data)][venue] = 1
							else:
								authors_venues_years[year_key][str(auth.firstChild.data)][venue] += 1
					if len(node.getElementsByTagName('author')) > 1:
						authors = [
							remap_name(a.firstChild.data)
							for a in node.getElementsByTagName('author')
						]
					if len(authors) > 1:
						authors.sort()
					for a, b in itertools.combinations(authors, 2):
						src = min(a,b)
						dst = max(a,b)
						if src not in existing_edges_years_alg[year_key].keys():
							existing_edges_years_alg[year_key][src] = set([dst])
							existing_edges_years[year_key][src] = [dst]
							#print(f"{src} {dst}", file=output)
						elif dst not in existing_edges_years_alg[year_key][src]:
							existing_edges_years_alg[year_key][src].add(dst)
							existing_edges_years[year_key][src].append(dst)
							#print(f"{src} {dst}", file=output)
	#output.close()
	dict_mapping = {"authormapping": names, "authorsVenuesYears": authors_venues_years, "graphEdgesYear": existing_edges_years, "yearLegend": years, "yearMap":year_map}
	json_obj = json.dumps(dict_mapping)
	with open(local_mapping_json, "w") as outj:
		outj.write(json_obj)

if __name__ == '__main__':
	#dataset_dblp()
	dataset_temporal_dblp()
