from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QFileDialog, QButtonGroup, QRadioButton
import gzip
from Bio import SeqIO
from collections import Counter

#common = ATCAGTGACTTC

class MyWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle('NGS Data Analysis Code 07 September 2023')
        layout = QVBoxLayout()
        
        # Left section (Text inputs)
        left_section = QVBoxLayout()        
        # Text inputs
        self.text_inputs = [
            #{'label': 'How many amino acids is the target sequence?', 'input': QLineEdit()},
            {'label': 'What is the conserved sequence downstream of the target codons?','input': QLineEdit()},
            #{'label': 'What is the 5 amino acid sequence that you are interested in tracking?', 'input': QLineEdit()},
           # {'label': 'What position is that amino acid in?','input':QLineEdit()},
            {'label': 'What is the name of your output file?','input':QLineEdit()},

        ]

        for text_input in self.text_inputs:
            left_section.addWidget(QLabel(text_input['label']))
            left_section.addWidget(text_input['input'])
        
        layout.addLayout(left_section)
        
        # Right section (File browsers)
        right_section = QVBoxLayout()

        # Labels, inputs, and buttons for file browsers
        self.file_browsers = [
            {'label': 'FASTQ.gz Data File B1:', 'input': QLineEdit(), 'button': QPushButton('Browse for B1 data')},
            {'label': 'FASTQ.gz Data File B2:', 'input': QLineEdit(), 'button': QPushButton('Browse for B2 data')},
            {'label': 'FASTQ.gz Data File B3:', 'input': QLineEdit(), 'button': QPushButton('Browse for B3 data')},
            {'label': 'FASTQ.gz Data File B4:', 'input': QLineEdit(), 'button': QPushButton('Browse for B4 data')}
        ]

        for file_browser in self.file_browsers:
            file_browser['button'].clicked.connect(self.browse_file)
            right_section.addWidget(QLabel(file_browser['label']))
            right_section.addWidget(file_browser['input'])
            right_section.addWidget(file_browser['button'])
            
        layout.addLayout(right_section)

        # Button to submit inputs
        self.button_submit = QPushButton('Submit')
        self.button_submit.clicked.connect(self.submit_inputs)
        layout.addWidget(self.button_submit)

        self.setLayout(layout)
        self.show()
        
        
        

    def browse_file(self):
        file_dialog = QFileDialog()
        if file_dialog.exec_():
            file_path = file_dialog.selectedFiles()[0]
            sender = self.sender()
            for file_browser in self.file_browsers:
                if file_browser['button'] is sender:
                    file_browser['input'].setText(file_path)
                    break

    def submit_inputs(self):
        self.text_input_values = [text_input['input'].text() for text_input in self.text_inputs]
        self.file_input_values = [file_browser['input'].text() for file_browser in self.file_browsers]
        
        fasta_files = [file_browser['input'].text() for file_browser in self.file_browsers]

        output_file = self.text_input_values[1]  # Assuming output file name is the second text input

        common = self.text_input_values[0]  # Assuming common sequence is the first text input
    
        # Step 2: Perform analysis on the Fasta file
        result_list = []
        for num, fasta_file in enumerate(fasta_files):
            [common_count, num_reads, AA_countv2] = fasta_to_deg_dict(fasta_file, common)
        #Create the result list
            for key, value in AA_countv2.items():
                found = False
                for entry in result_list:
                    if entry[0] == key:
                        found = True
                        entry[num+1] = value   
                        break
             
                if not found:
                    inner_list = [key] + [0] * 5
                    inner_list[num+1] = value
                    result_list.append(inner_list)
        #Calculate PSI
            new_result_list = PSI_calc(result_list)
        write_output(output_file + f"_B{num + 1}.txt", common_count, num_reads, AA_countv2, new_result_list)            
        print( "Output Generated")
        
    
def user_inputs():
    FASTQ_file = window.file_input_values[0]
    output_file = str(window.text_input_values[1])+'.txt'
    common = window.text_input_values[0]
    #num_degen = int(window.text_input_values[0])
    #AA_target = window.text_input_values[2]
    #AA_position = int(window.text_input_values[3])
    return [FASTQ_file, output_file, common]

def fastq_to_fasta(input_file):
    input_split = input_file.split('.', 1)
    output_file = input_split[0] + '.fasta'
    with gzip.open(input_file, 'rt') as gzipped:
        # Open the output FASTA file
        with open(output_file, 'w') as fasta_out:
            # Parse each record from the FASTQ.gz file and write it to the FASTA file
            for record in SeqIO.parse(gzipped, 'fastq'):
                seq = str(record.seq).replace('\n', '')
                fasta_out.write(f'>{record.id}\n{seq}\n')
    return output_file

#Genetic code library
codon_table = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S', 'TCA':'S','TCG':'S','TAT':'Y','TAC':'Y', 'TAA':'*','TAG':'*','TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
                'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'CCT':'P', 'CCC':'P', 'CCA':'P','CCG':'P','CAT':'H','CAC':'H', 'CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
                'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
                'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}

def fasta_to_deg_dict(FASTA_files, common):
    results = []
    common_count = 0
    num_reads = 0
    temp_codons_list = []
    AA_list = []    
    for line in open(FASTA_files, 'r'):
            if line[0] != ">":
                num_reads += 1
                common_index = line.find(common)
                if common_index >=15:
                    common_count += 1
                    temp_codons_list.append(line[common_index-15:common_index])
                

    for entry in temp_codons_list: 
        protein = ''
        for i in range (0, 15, 3):
            codon = entry[i:i+3]
            amino_acid = codon_table.get(codon, "X")
            protein += amino_acid
        AA_list.append(protein)
    AA_count = Counter(AA_list)
    AA_countv2 = dict(sorted(AA_count.items(), key=lambda item: item[1], reverse=True)) 
    results.append([common_count, num_reads, AA_countv2])
    return [common_count, num_reads, AA_countv2]

def PSI_calc(result_list, RC, Nt_AA):
    for entry in result_list:
        totalreads = entry[1] + entry[2] + entry[3] + entry[4]
        if totalreads >= int(RC) and entry[0][0] == str(Nt_AA):
            PSI = (4*entry[1]/totalreads)+(3*entry[2]/totalreads)+(2*entry[3]/totalreads)+(1*entry[4]/totalreads)
            entry[5] = round(PSI, 2)
    new_result_list = sorted(result_list, key=lambda x: x[5], reverse=True)
    return new_result_list


def write_output(output_file, common_count, num_reads, AA_count, result_list):
    '''/
    Output information to file
    '''
    #AA_count_common = AA_count.most_common()
 
    with open(output_file,'w') as output_file:
        output_file.write('Across 4 bins of data, here are the read counts and PSI for each 5 AA sequence.\n\n')
        #output_file.write('The number of times the common, post-degenerate codon sequence appears is '+ str(common_count)+ ' times across '+ str(int(num_reads-1))+ ' reads.\n\n')
        #output_file.write('You are interested in '+str(AA_target)+ ' at position '+str(AA_position)+'. It appears '+str(target_count)+ ' times in this file.\n\n')
        #for key, value in AA_count.items():
            #output_file.write(str(key) + ',' + str(value) + '\n')
        for entry in result_list:
            if '*' not in entry[0] and (entry[1]+entry[2]+entry[3]+entry[4]) >24 and entry[0][0] == 'F':
                output_file.write('Sequence: ' + str(entry[0]) + ', B1: ' + str(entry[1]) + ', B2:' + str(entry[2])+ ', B3:' + str(entry[3]) + ', B4:' + str(entry[4]) + ' PSI: ' +str(entry[5]) +'\n')
                #output_file.write((str(entry[0])+'\n')*(entry[1]+entry[2]+entry[3]+entry[4]))


if __name__ == '__main__':
    app = QApplication([])
    window = MyWindow()
    app.exec_()
    [NGS_files, output_file, common] = user_inputs()
    write_output(output_file + f"_output_B{idx}.txt", common_count, num_reads, AA_countv2, result_list)
    print(f"B{idx} Output Generated")
    print ("Run completed.")
