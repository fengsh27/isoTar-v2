import numpy as np
import sys, os, argparse, getopt, pdb
from itertools import islice

# -----------------------------------------------
def writeFile(data, filename, mode="w"):
    d = os.path.dirname(filename)

    if d != "":
        if not os.path.exists(d):
            os.makedirs(d)

    fl = open(filename, mode)
    fl.write(data)
    fl.close()

# -----------------------------------------------
def writeDataTableAsText(data, filename, mode="w"):
    text = formatDataTable(data, "\t", "\n")

    writeFile(text, filename, mode)

# -----------------------------------------------
def formatDataTable(data, col_sep="\t", row_sep="\n"):
    return row_sep.join([col_sep.join([str(item1) for item1 in item]) for item in data])

# -----------------------------------------------------------------
def addPadding(seq, length):

    total_pad = abs(len(seq)-length)
    left_pad = total_pad//2
    right_pad = total_pad-left_pad

    if len(seq) > length:
        return seq[left_pad:left_pad+length]

    return seq + ''.join(['N']*total_pad)
    #return ''.join(['N']*left_pad) + seq + ''.join(['N']*right_pad)

# -----------------------------------------------------------------
def oneHotEncode(seq):

    nt_map = {'A':0, 'T':1, 'U':1, 'C':2, 'G':3}

    code = np.zeros((len(seq), 4))

    for i in range(len(seq)):
        if seq[i] == 'N': code[i] = [0.25]*4
        else:
            assert seq[i] in nt_map, '{} not recognized for one hot coding'.format(seq[i])
            code[i, nt_map[seq[i]]] = 1

    return code

#------------------------------------------------------------------------------
def processData(interactions):

    X_mirnas, X_targets = [], []

    for mirna_id, target_id, mirna_seq, target_seq in interactions:

        mirna_seq_ext = addPadding(mirna_seq, len_mirna)
        mirna_code = oneHotEncode(mirna_seq_ext)
        target_seq_ext = addPadding(target_seq, len_target)
        target_code = oneHotEncode(target_seq_ext)
        X_mirnas.append(mirna_code)
        X_targets.append(target_code)

    return X_mirnas, X_targets

#------------------------------------------------------------------------------
def loadModel(model_path):

    json_file = open(model_path + '.json', 'r')
    loaded_model_json = json_file.read()
    json_file.close()
    model = model_from_json(loaded_model_json)
    # load weights into new model
    model.load_weights(model_path + ".h5")
    print("Loaded model from disk")
    return model

# -------------------------------------------------------------------------
def showArgError(parser):
    parser.parse_args(['-h'])
    exit()

# -------------------------------------------------------------------------
def showError(msg):
    print('----------- Error !!! ------------')
    print(msg)
    print('----------------------------------')

# -------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Predicts interactions between miRNA and target pair. ' + \
                                             'It takes an input file of query miRNA/isomiR and target pairs and predict the interactions between the pairs. ' + \
                                             'Alternatively, it can take separate files for miRNA/isomiR and target sequences and outputs all the possible interactions.', \
                                 epilog='Example: "python3 dmiso.py -p examples/test_pairs.txt -o examples/test_output.txt" OR ' + \
                                                  '"python3 dmiso.py -m examples/test_miRNAs.fa -t examples/test_mRNAs.fa -o examples/test_output.txt"')
required_arg = parser.add_argument_group('required arguments')
required_arg.add_argument('-p', help="Path for miRNA (or isomiR) and mRNA pair file in a tsv format. The file must contain data in the following order: miRNA id, target id, miRNA sequence, target sequence", required=True,
                          metavar="PAIR")
required_arg.add_argument('-m', help="Path for miRNA or isomiR sequence in fasta format.", required=True, metavar="MIRNA/ISOMIR")
required_arg.add_argument('-t', help="Path for mRNA sequence in fasta format.", required=True, metavar="MRNA")
parser.add_argument('-o', help='Path for DMISO outputs', metavar="OUTPUT")

if '-h' in sys.argv[1:]:
    showArgError(parser)

opts, args = getopt.getopt(sys.argv[1:], 'p:m:t:o:')

pair_file_path = ''
mi_seq_file_path = ''
m_seq_file_path = ''
output_path = ''

if len(opts) == 0: showArgError(parser)

for i in opts:
    if i[0] == '-p':
        pair_file_path = i[1]
        if pair_file_path == '': showArgError(parser)
    elif i[0] == '-m':
        mi_seq_file_path = i[1]
        if mi_seq_file_path == '': showArgError(parser)
    elif i[0] == '-t':
        m_seq_file_path = i[1]
        if m_seq_file_path == '': showArgError(parser)
    elif i[0] == '-o':
        output_path = i[1]
        if output_path == '': showArgError(parser)
    else:
        showArgError(parser)

is_valid_pair_path = os.path.isfile(pair_file_path)
is_valid_mi_seq_path = os.path.isfile(mi_seq_file_path)
is_valid_m_seq_path = os.path.isfile(m_seq_file_path)

# Validate pair path
if not (is_valid_pair_path or (is_valid_mi_seq_path and is_valid_m_seq_path)):
    if not is_valid_pair_path:
        showError('The following file path for miRNA target pair file does not exist' + '\n' + pair_file_path)
    elif not is_valid_mi_seq_path:
        showError('The following file path for miRNA (isomiR) sequence file does not exist' + '\n' + mi_seq_file_path)
    elif not is_valid_m_seq_path:
        showError('The following file path for mRNA sequence file does not exist' + '\n' + m_seq_file_path)
    exit()

if output_path == '': output_path = os.path.join(os.path.dirname(pair_file_path), 'predictions_dmiso.txt')

from keras.models import model_from_json

tool_path = os.path.dirname(os.path.realpath(__file__))
model_path = os.path.join(tool_path, 'models/model_top')
model = loadModel(model_path)
_, len_mirna, _ = model.get_layer('mirna').output_shape
_, len_target, _ = model.get_layer('target').output_shape

header = [['miRNA ID', 'Target ID', 'miRNA Sequence', 'Target Sequence', 'Prediction Score', 'Prediction']]
writeDataTableAsText(header, output_path)

if is_valid_pair_path:

    with open(pair_file_path) as f:
        while True:
            file_chunk = list(islice(f, 1000000))
            if not file_chunk:
                break

            interactions = []
            for item in file_chunk:
                item = item.strip()
                mi_id, m_id, mi_seq, m_seq = item.split('\t')
                while len(m_seq) > len_target * 1.5:
                    m_sub_seq = m_seq[:len_target]
                    m_seq = m_seq[len_target:]
                    interactions.append([mi_id, m_id, mi_seq, m_sub_seq])
                interactions.append([mi_id, m_id, mi_seq, m_seq])

            X_mirnas, X_targets = processData(interactions)
            y_pred_pp = model.predict([X_mirnas, X_targets])[:, 0]
            y_pred = (y_pred_pp > 0.5).astype(int)
            output = np.c_[interactions, y_pred_pp, y_pred]

            writeFile("\n", output_path, "a")
            writeDataTableAsText(output, output_path, "a")

elif is_valid_mi_seq_path and is_valid_m_seq_path:

    chunk_size = 2000000

    file_chunk_mi = []

    with open(mi_seq_file_path) as f_mi:
        while True:
            file_chunk_mi += list(islice(f_mi, chunk_size))
            if not file_chunk_mi:
                break

            if len(file_chunk_mi) >= chunk_size:
                for i in range(len(file_chunk_mi) - 1, -1, -1):
                    if file_chunk_mi[i][0] == '>':
                        break
                seq_info = ''.join(file_chunk_mi[:i]).split('>')
                file_chunk_mi = file_chunk_mi[i:]
            else:
                seq_info = ''.join(file_chunk_mi).split('>')
                file_chunk_mi = []

            mi_info = {}
            for item in seq_info[1:]:
                lines = item.split('\n')
                id = lines[0]
                seq = ''.join(lines[1:-1])
                mi_info[id] = seq.upper()

            file_chunk_m = []

            with open(m_seq_file_path) as f_m:
                while True:
                    file_chunk_m += list(islice(f_m, chunk_size))
                    if not file_chunk_m:
                        break

                    if len(file_chunk_m) >= chunk_size:
                        for i in range(len(file_chunk_m) - 1, -1, -1):
                            if file_chunk_m[i][0] == '>':
                                break
                        seq_info = ''.join(file_chunk_m[:i]).split('>')
                        file_chunk_m = file_chunk_m[i:]
                    else:
                        seq_info = ''.join(file_chunk_m).split('>')
                        file_chunk_m = []

                    m_info = {}
                    for item in seq_info[1:]:
                        lines = item.split('\n')
                        id = lines[0]
                        seq = ''.join(lines[1:-1])
                        m_info[id] = seq.upper()

            interactions = []
            for mi_id, mi_seq in mi_info.items():
                for m_id, m_seq in m_info.items():
                    while len(m_seq) > len_target*1.5:
                        m_sub_seq = m_seq[:len_target]
                        m_seq = m_seq[len_target:]
                        interactions.append([mi_id, m_id, mi_seq, m_sub_seq])
                    interactions.append([mi_id, m_id, mi_seq, m_seq])
            X_mirnas, X_targets = processData(interactions)
            y_pred_pp = model.predict([X_mirnas, X_targets])[:, 0]
            y_pred = (y_pred_pp > 0.5).astype(int)
            output = np.c_[interactions, y_pred_pp, y_pred]
            output = output[y_pred == 1]

            writeFile("\n", output_path, "a")
            writeDataTableAsText(output, output_path, "a")

print('Prediction complete !!')
print('Please check the prediction results at ' + os.path.abspath(output_path))