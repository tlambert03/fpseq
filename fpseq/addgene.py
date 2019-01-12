import os
import io
import requests
import math
import warnings
import pandas as pd
import sqlite3
import re
import json
from subprocess import run, PIPE
from bs4 import BeautifulSoup
from Bio import SeqIO, BiopythonParserWarning
from collections import defaultdict, Counter, OrderedDict
from fpseq.mutations import get_mutations

warnings.simplefilter('ignore', BiopythonParserWarning)

# how many labeled m have the 206K mutation


FILESTORE = os.path.expanduser('/Users/talley/Dropbox (HMS)/addgene')
#db = sqlite3.connect(os.path.join(FILESTORE, 'database.db'))
#cursor = db.cursor()
# TABLE: sequences
#   sequence_id
#   label
#   date
#   contributor
#   length
#   cds
#   plasmid_id

# if not cursor.execute("SELECT name FROM sqlite_master WHERE type='table';").fetchall():
#     db.execute("""CREATE TABLE plasmids
#                (plasmid_id integer PRIMARY KEY, name text, type text, depositor text, tags text,
#                 expression text, mutation text, flame text, purpose text, use text);""")
#     db.execute(""" CREATE TABLE sequences (
#                sequence_id integer,
#                label text,
#                date text,
#                contributor text,
#                length text,
#                cds text,
#                plasmid_id integer NOT NULL,
#                         FOREIGN KEY (plasmid_id) REFERENCES plasmids(plasmid_id),
#                         UNIQUE(sequence_id, label)
#                );""")


def get_addgene_ids(query, withtags=True):
    if withtags:
        url = 'https://www.addgene.org/search/advanced/?q={0}&tags={0}&results_per_page=20&page={1}&sort_type=popularity'
    else:
        url = 'https://www.addgene.org/search/advanced/?q={0}&results_per_page=20&page={1}&sort_type=popularity'
    response = requests.get(url.format(query, 1))
    soup = BeautifulSoup(response.text, 'lxml')
    count = soup.find('span', {'id': 'results_count'})
    try:
        pages = math.ceil(int(count.text.split(' of ')[1]) / 20)
    except Exception:
        pages = 1

    ids = {}
    for page in range(1, pages + 1):
        print("page: {} of {}".format(page, pages))
        if page > 1:
                response = requests.get(url.format(query, page))
                soup = BeautifulSoup(response.text, 'lxml')
        for result in soup.find_all('div', {'class': 'search-result-plasmid'}):
            # (id, name, type, depositor, tags)
            container = result.find('div', {'class': 'search-results-container'})
            itemtype = result.find('span', {'class': 'addgene-item-type'})

            plasID = result.attrs['id'].replace('search_result_', '')
            name = result.find('div', {'class': 'search-result-plasmid-title'}).find('a').text
            ids[plasID] = {
                'name': name,
                'type': itemtype.text if itemtype else None,
                'depositor': None,
                'expression': None,
                'tags': None,
                'mutation': None,
                'flame': None,
                'purpose': None,
                'use': None
            }

            for flamelev in ('low', 'medium', 'high'):
                if result.find('span', {'class': 'addgene-flame-{}'.format(flamelev)}):
                    ids[plasID]['flame'] = flamelev

            rows = container.find_all('div', {'class': 'row'})

            for row in rows:
                try:
                    lab1 = row.find_all('div')[0].text.strip()
                    if lab1 in ('Expression', 'Use', 'Purpose'):
                        ids[plasID][lab1.lower()] = row.find_all('div')[1].text.strip()
                except Exception:
                    pass
                try:
                    lab2 = row.find_all('div')[2].text.strip()
                    if lab2 in ('Depositor', 'Tags', 'Mutation'):
                        ids[plasID][lab2.lower()] = row.find_all('div')[3].text.strip()
                except Exception:
                    pass

    return ids


def get_addgene_genbanks(ids, force=False, filestore=FILESTORE):

    for id in ids:  # (id, name, type, depositor, tags)
        if any(['plasmid-{}'.format(id) in x for x in os.listdir(filestore)]):
            if not force:
                print("id {:6} already downloaded... skipping".format(id))
                continue
        else:
            print(id, " not found")
        response = requests.get('https://www.addgene.org/{}/sequences/'.format(id))
        soup = BeautifulSoup(response.text, 'lxml')

        seqIDs = []
        for dep in ('depositor', 'addgene'):
            seqtype = '{}-full'.format(dep)
            section = soup.find('section', {'id': seqtype})
            if not section:
                seqtype = '{}-partial'.format(dep)
                section = soup.find('section', {'id': seqtype})
            if section:
                link = section.find('a', {'class': 'genbank-file-download'}).attrs['href']
            else:
                continue

            filename = link.split('/')[-1].replace('.gbk', '-{}.gbk'.format(seqtype))
            seqIDs.append(filename.split('sequence-')[1].split('-')[0])

            os.mkdir(filestore) if not os.path.exists(filestore) else None

            if filename in os.listdir(filestore):
                print("found {}".format(filename))
                with open(os.path.join(filestore, filename), 'r') as f:
                    content = f.read()

            if filename not in os.listdir(filestore):
                content = requests.get(link)
                content = content.text if content else None
                with open(os.path.join(filestore, filename), 'w') as f:
                    f.write(content)
                print("wrote {}".format(filename))


def list_CDS_variants(parsed_records):
    out = defaultdict(lambda: defaultdict(set))
    for key, value in parsed_records.items():
        for feature in value['features']:
            name = feature.qualifiers.get('label', [''])[0].strip('')
            trans = feature.qualifiers.get('translation', [''])[0].strip('')
            print(feature.qualifiers)
            if name and trans:
                out[name][trans].add(key)
    return {k: dict(v) for k, v in out.items()}


def align_list(seq_list, binary='muscle', output='clw'):

    def fastas(seqs):
        return io.StringIO("\n".join([">SEQ_{}\n{}".format(n, s) for n, s in enumerate(seqs)]))

    fasta = fastas(seq_list)
    cmd = [binary]
    # faster
    cmd += ['-maxiters', '6', '-diags', '-quiet', '-' + output]
    result = run(cmd, input=fasta.read(), stdout=PIPE, encoding='ascii')
    return result.stdout


def slugify(value, allow_unicode=False):
    value = str(value)
    value = re.sub(r'[^\w\s-]', '', value).strip().lower()
    return re.sub(r'[-\s]+', '-', value)


def from_fpbase(slug):
    url = 'https://www.fpbase.org/api/{}/?format=json'.format(slugify(slug))
    response = requests.get(url)
    return json.loads(response.content).get('seq')


class AddgeneDB(object):
    def __init__(self, filestore=None, dbname='database'):
        if filestore is None:
            filestore = os.path.expanduser('/Users/talley/Dropbox (HMS)/addgene')
        self.filestore = filestore
        self.db = sqlite3.connect(os.path.join(filestore, dbname + '.db'))
        self.cursor = self.db.cursor()
        self._fpbase_seqs = {}

    def commit(self):
        self.db.commit()

    def close(self):
        self.db.close()

    def sql(self, query):
        return pd.read_sql_query(query, self.db)

    def all_data(self):
        return pd.read_sql_query('select * from sequences inner join plasmids on '
                                 'plasmids.plasmid_id = sequences.plasmid_id', self.db)

    def most_common_variants(self, name, depth=1, with_counts=False):
        variants = self.get_variants(name)
        if variants:
            if with_counts:
                out = [(i[0], len(i[1])) for i in sorted(variants.items(),
                       key=lambda kv: len(kv[1]), reverse=True)[:depth]]
            else:
                out = [i[0] for i in sorted(variants.items(),
                       key=lambda kv: len(kv[1]), reverse=True)[:depth]]
            return out

    def get_addgene_ids(self, query, withtags=True):
        ids = get_addgene_ids(query, withtags=withtags)
        self.cursor.executemany('INSERT OR REPLACE INTO plasmids VALUES (?,?,?,?,?,?,?,?,?,?)',
                                [(int(k), d['name'], d['type'], d['depositor'], d['tags'],
                                 d['expression'], d['mutation'], d['flame'], d['purpose'], d['use'])
                                 for k, d in ids.items()])
        self.db.commit()
        return ids

    def to_csv(self):
        self.cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = self.cursor.fetchall()
        for table_name in tables:
            table_name = table_name[0]
            table = pd.read_sql_query("SELECT * from %s" % table_name, self.db)
            table.to_csv(os.path.join(
                self.filestore, 'db_' + table_name + '.csv'), index_label='index')

    def counts(self):
        T = self.cursor.execute('SELECT label, cds FROM sequences').fetchall()
        A = set(T)
        c = Counter([a[0] for a in A])
        t = Counter([t[0] for t in T])
        for name, count in c.most_common():
            print("{:3} variants in {:4} total: {}".format(count, t[name], name))
        print()
        print('Most Variable')
        ratios = [(name, count / t[name]) for name, count in c.items() if (t[name] > 1 and count > 1)]
        ratios.sort(key=lambda x: x[1], reverse=True)
        for name, ratio in ratios:
            print("{:4}% - {}".format(round(ratio * 100), name))
        # return c

    def search_plasmids_for_sequence(self, seq):
        if len(seq) < 10:
            seq = self._get_fpbase_seq(seq)
        d = self.cursor.execute('SELECT plasmid_id FROM sequences WHERE cds LIKE "%{}%"'.format(seq)).fetchall()
        return list(set([x[0] for x in d]))

    def get_variants(self, label, exact=True):
        '''Return a dict of all the variants (and which plasmid_ids have them)
        for a given label'''
        if exact:
            r = self.cursor.execute('SELECT cds, plasmid_id FROM sequences WHERE label = ?', [label]).fetchall()
        else:
            r = self.cursor.execute('SELECT cds, plasmid_id FROM sequences WHERE label LIKE "%{}%"'.format(label)).fetchall()
        D = defaultdict(set)
        for x in r:
            D[x[0]].add(x[1])
        return D

    def variant_mutations(self, label, reference='avGFP', ignore_ends=True, exact=True):
        """ find sequences for a given label that vary from a reference sequence """
        v = self.get_variants(label, exact=exact)
        if len(reference) < 10:
            reference = self._get_fpbase_seq(reference)
        muts = defaultdict(set)
        for variant, ids in v.items():
            ms = get_mutations(reference, variant)
            if ignore_ends:
                removes = [x for x in ms if (x.start_idx < 4 or x.start_idx > 233)]
                for r in removes:
                    ms.muts.remove(r)
            muts[ms].update(ids) if ms else None
        return OrderedDict(sorted(muts.items(), key=lambda x: len(x[1]), reverse=True))

    def _get_fpbase_seq(self, seq):
        seq = slugify(seq)
        if seq in self._fpbase_seqs:
            return self._fpbase_seqs[seq]
        else:
            r = from_fpbase(seq)
            self._fpbase_seqs[seq] = r
            return r

    def get_mutations(self, mutation='A206K', ref='avGFP', names=None, without=False):
        """ searches the database for sequences that possess (or lack) a certain mutation
        if without == True, returns hits that do NOT have the mutation
        """
        ref = self._get_fpbase_seq(ref)

        if names == 'all':
            data = self.cursor.execute('SELECT label, cds, plasmid_id FROM sequences').fetchall()
        elif names is not None:
            if isinstance(names, str):
                names = [names]
            data = self.cursor.execute('SELECT label, cds, plasmid_id FROM sequences'
                                       ' WHERE label IN (%s)' % ','.join('?' * len(names)), names).fetchall()
        else:  # get names starting with m
            data = self.cursor.execute('SELECT label, cds, plasmid_id FROM sequences'
                                       ' WHERE label like "m%"').fetchall()

        hits = defaultdict(set)
        for label, cds, plasmid_id in data:
            mutset = get_mutations(ref, cds)
            if without:
                if mutation not in mutset:
                    hits[label].add(plasmid_id)
            else:
                if mutation in mutset:
                    hits[label].add(plasmid_id)

        for k, v in hits.items():
            q = self.cursor.execute('SELECT tags, plasmid_id from plasmids where plasmid_id'
                                    ' in (%s)' % ','.join('?' * len(v)), [str(x) for x in v]).fetchall()
            for tags, plasmid_id in q:
                if tags and (k in tags):
                    print("{}: {}".format(k, plasmid_id))
        return hits

    def find_sequence_inconsistencies(self):
        self.cursor.execute('select cds, contributor, label, length, plasmid_id from sequences')
        data = self.cursor.fetchall()
        d = defaultdict(list)
        for cds, contributor, label, length, plasmid_id in data:
            d[plasmid_id].append((cds, contributor, label, length, plasmid_id))
        for k, v in d.items():
            if len(v) > 1:  # if it has multiple coding sequences
                uniques = set([(x[0], x[2]) for x in v])
                counter = Counter([u[1] for u in uniques])
                if counter.most_common()[0][1] > 1:
                    lab = counter.most_common()[0][0]
                    seqs = [(i[0], i[1], i[3]) for i in v if i[2] == lab]
                    ms = get_mutations(seqs[0][0], seqs[1][0])
                    removes = [x for x in ms if (x.start_idx < 3
                               or x.start_idx > len(seqs[0][0]) - 3)]
                    for r in removes:
                        ms.muts.remove(r)
                    if ms:
                        print("{} in plasmid {} - {} {} -> {} {} - mutations: {}"
                              .format(lab, k, seqs[0][1], seqs[0][2],
                                      seqs[1][1], seqs[1][2], ms))

    def get_addgene_genbanks(self, ids=None):
        if ids is None:
            L = self.cursor.execute('SELECT plasmid_id FROM plasmids').fetchall()
            ids = [x[0] for x in L]
        get_addgene_genbanks(ids, filestore=self.filestore)

    def parse_genbanks(self, ids=None, minlength=140):
        ignoredCDS = ('NeoR', 'AmpR', 'KanR', '6xHis', 'Xpress(TM)', 'T7 tag',
                      'HygR', 'PuroR', 'Rep101', 'araC', 'GmR', 'TcR', 'TetR',
                      'ccdB', 'BSD', 'pVS1', 'CmR', 'URA3', 'ORF1629', 'ORF603',
                      'SmR', 'VN155', 'LEU2', 'Cre', 'NrsR', 'LAMP1', 'DHFR',
                      'BirA', 'RhoA', 'mini-white', 'rtTA', 'lacZ', 'c-Myc', 'TRP1',
                      'URA3', 'HIS3', 'FLP', 'Cas9', 'dCas9', 'lacI', 'VN173',
                      't antigen', 'BlpR')
        filelist = [f for f in os.listdir(self.filestore) if f.endswith('.gbk')]
        if isinstance(ids, (list, dict, tuple, set)):
            filelist = [f for f in filelist if any([str(i) in f for i in ids])]

        for file in filelist:
            with open(os.path.join(self.filestore, file), 'r') as gbfile:
                try:
                    record = list(SeqIO.parse(gbfile, "genbank"))[0]
                except Exception:
                    print(file, ' failed')
            seqID = int(file.split('sequence-')[1].split('-')[0])
            plasID = int(file.split('plasmid-')[1].split('-')[0])
            date = record.annotations.get('date')
            contrib = 'depositor' if 'depositor' in file else 'addgene'
            length = 'full' if 'full' in file else 'partial'

            cds = []
            for sf in record.features:
                if (sf.type == 'CDS'
                        and (len(sf.qualifiers.get('translation', [[]])[0]) > minlength)
                        and not any([x in sf.qualifiers.get('label', [None])[0] for x in ignoredCDS])):
                    cds.append(sf)
            cds.sort(key=lambda x: len(x.qualifiers.get('translation', [[]])[0]), reverse=True)
            for cd in cds:
                trans = cd.qualifiers.get('translation', [[]])[0]
                label = cd.qualifiers.get('label', [None])[0]

                self.cursor.execute('INSERT OR REPLACE INTO sequences VALUES (?,?,?,?,?,?,?)',
                                    (seqID, label, date, contrib, length, trans, plasID))
                self.db.commit()


def addgene_check():
    from proteins.models.lineage import Lineage
    db = AddgeneDB()

    good = []
    missing = []
    print('protein             muts to addgene')
    print('-----------------------------------')
    for L in Lineage.objects.all().select_related('protein'):
        v = db.most_common_variants(L.protein.name)
        if not v:
            missing.append(L.protein.name)
            continue
        if L.protein.seq == v[0]:
            good.append(L.protein.name)
        if L.protein.seq != v[0]:
            muts = L.protein.seq.mutations_to(v[0])
            if muts:
                print(f'{L.protein.name:<20}{muts}')  # noqa
            else:
                print(L.protein.name)
                print(f'FPb: {L.protein.seq}\nAdg: {v[0]}')
    print("\nGOOD")
    print(", ".join(good))
    print("\nMISSING")
    print(",".join(missing))

# mdb = AddgeneDB()
# Qs = set(('Keima', 'Tomato', 'Wasabi', 'Emerald', 'Maple', 'Dendra', 'Cerulean', 'Ruby',
#           'Scarlet', 'TFP', 'Plum', 'Papaya', 'Orange', 'Neptune', 'Grape', 'Nectarine',
#           'mKO', 'Cyan', 'Keima', 'Kusabira', 'Midoriishi', 'IFP', 'Honeydew', 'DsRed',
#           'Dronpa', 'Clover', 'AsRed', 'Citrine', 'Dreiklang', 'CyPet', 'Ypet', 'Katushka',
#           'Eos', 'Turquoise', 'mStable', 'mGeos', 'Garnet', 'Blueberry', 'Azami', 'mClav',
#           'Carmine', 'Cardinal', 'Kate', 'Apple', 'Strawberry', 'Crimson', 'Banana',
#           'eqFP', 'amilFP', 'amilCP', 'Sapphire', 'Skylan', 'CaMP', 'KillerRed',
#           'Kaede', 'HcRed', 'Gamillus', 'FusionRed', 'Clomeleon', 'Aquamarine',
#           'Azurite', 'Kalama', 'Lumin', 'Maroon', 'Raspberry', 'Tangerine', 'Topaz',
#           'Ultramarine', 'amFP', 'ZsGreen', 'SuperNova', 'Sandercyanin', 'SHardonnay',
#           'Padron', 'LanFP', 'KikG', 'Kohinoor', 'CyOFP', 'Ametrine', 'NeonGreen',
#           'IrisFP', 'mKG', 'iFP1', 'iFP2', 'Sirius', 'mKO', 'GFP', 'RFP', 'YFP'))
# ids = {}
# for q in Qs:
#     print(q)
#     ids.update(mdb.get_addgene_ids(q, withtags=False))
# get_addgene_genbanks(ids)
# parse_genbanks(ids)
