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
        self.setWindowTitle('NGS Data Analysis Code 17 May 2023')
        layout = QVBoxLayout()
        
        # Left section (Text inputs)
        left_section = QVBoxLayout()        
        
        # Text inputs
        self.text_inputs = [
            {'label': 'How many amino acids is the target sequence?', 'input': QLineEdit()},
            {'label': 'What is the conserved sequence downstream of the target codons?', 'input': QLineEdit()},
            {'label': 'What is the first amino acid you are interested in tracking?', 'input': QLineEdit()},
            {'label': 'What position is the first amino acid in?', 'input': QLineEdit()},
            {'label': 'What is the second amino acid you are interested in tracking?', 'input': QLineEdit()},
            {'label': 'What position is the second amino acid in?', 'input': QLineEdit()},
            {'label': 'What is the name of your output file?', 'input': QLineEdit()},
        ]

        for text_input in self.text_inputs:
            left_section.addWidget(QLabel(text_input['label']))
            left_section.addWidget(text_input['input'])
        
        layout.addLayout(left_section)
        
        # Right section (File browsers)
        right_section = QVBoxLayout()

        # Labels, inputs, and buttons for file browsers
        self.file_browsers = [
            {'label': 'FASTQ.gz Data File:', 'input': QLineEdit(), 'button': QPushButton('Browse for B1 data')}
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

        [NGS_file, common, num_degen, AA_target1, AA_position1, AA_target2, AA_position2, output_file] = user_inputs()
    
        # Step 1: Convert FastQ.gz to Fasta
        fasta_file = fastq_to_fasta(NGS_file)
        print("FASTA Generated")
    
        # Step 2: Perform analysis on the Fasta file
        [degenlist, codons_list, common_count, num_reads] = get_degen_list(fasta_file, common, int(num_degen))
        Count_full_AA = translate_degen_list(degenlist, num_degen)
        target_count1 = position_counting1 (Count_full_AA, AA_target1, AA_position1, num_degen)
        target_count2 = position_counting2 (Count_full_AA, AA_target2, AA_position2, num_degen)
        [Positions, AA_counts], target_counter_both = before_after (AA_target1, AA_position1,AA_target2,AA_position2, num_degen, Count_full_AA )
        write_output(output_file, AA_target1, AA_position1,AA_target2, AA_position2, common_count, target_count1,target_count2,target_counter_both, num_reads, num_degen, AA_counts, Positions)
        print("Output Generated")
        
    
def user_inputs():
    FASTQ_file = window.file_input_values[0]
    output_file = str(window.text_input_values[6]) + '.txt'
    num_degen = int(window.text_input_values[0])  # Assuming num_degen is provided as is, without conversion to int
    common = window.text_input_values[1]
    AA_target1 = window.text_input_values[2]
    AA_position = int(window.text_input_values[3])
    AA_target2 = window.text_input_values[4]
    AA_position2 = int(window.text_input_values[5])
    return [FASTQ_file, common, num_degen, AA_target1, AA_position, AA_target2, AA_position2, output_file]

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


def get_degen_list(NGS_file, common, num_degen):
    '''
    Creates list of 'num_degen' degenerate codons that appear before the common sequence
    degenlist contains strings of all num_degen codons with no separation (e.g. ['TAGATGCATATTACT'])
    codons_list contains an array of each codon for a single line (e.g. [['TAG', 'ATG', 'CAT', 'ATT', 'ACT']])
    '''
    degenlist = []
    codons_list = []
    common_count = 0
    num_reads = 0
    for line in open(NGS_file, 'r'):
        if line[0] != ">":
            num_reads += 1
            common_index = line.find(common)
            #Find provides the position # in the string that the common sequence begins at (remember 0 is the first index)
            if common_index >= num_degen*3:
                #Excludes reads without a sufficient codon count
                common_count += 1
                temp_codons_list = []
                for i in reversed(range(1, num_degen+1)):
                    temp_codons_list.append(line[common_index-i*3:common_index-i*3+3])
                    #Adds [3i:3i+3] to temp_codons_list
                degenlist.append(line[common_index-num_degen*3:common_index])
                codons_list.append(temp_codons_list)
    return [degenlist, codons_list, common_count, num_reads]


#Generates a returned variable that is a dictionary, unsorted (i.e. [('TGA':'34')...] of amino acids
def translate_degen_list(degenlist,num_degen):
    AAlist = []
    for entry in degenlist:
        protein = ''
        for i in range (0, int(num_degen)*3, 3):
            codon = entry[i:i+3]
            amino_acid = codon_table.get(codon, "X")
            protein += amino_acid
        AAlist.append(protein)
    AA_count = Counter(AAlist)
    return AA_count
            
def position_counting1 (AA_count, AA_target1, AA_position1, num_degen):
    target_count1 = 0
    Pos2 = []
    if AA_position1 > num_degen:
        return []
    else:
        for key in AA_count:
            if key[int(AA_position1)-1] == str(AA_target1):
                target_count1 += int(AA_count[key])
    return target_count1      
      
def position_counting2 (AA_count, AA_target2, AA_position2, num_degen):
    target_count2 = 0
    Pos2 = []
    if AA_position2 > num_degen:
        return []
    else:
        for key in AA_count:
            if key[int(AA_position2)-1] == str(AA_target2):
                target_count2 += int(AA_count[key])
    return target_count2    


def before_after (AA_target1, AA_position1,AA_target2,AA_position2, num_degen, deg_library):
    target_count1 = 0
    target_count2 = 0
    target_counter_both = 0
    AA_match = ['A','C','D','E','F' ,'G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*']
    AA_counts = [{AA: 0 for AA in AA_match} for i in range(int(num_degen)+1)]

    for key in deg_library:
        if key[int(AA_position1)-1] == str(AA_target1) and key[int(AA_position2)-1] == str(AA_target2) :
        
            target_counter_both += int(deg_library[key])
            for i in range(len(key)):
                for AA in AA_match:
                    if i<int(AA_position1)-1 and key[i] == AA and i != int(AA_position2):
                        AA_counts[i][AA] += deg_library[key]
                    if i>int(AA_position1)-1 and key[i] == AA and i != int(AA_position2):
                        AA_counts[i-1][AA] += deg_library[key]
        
    Positions = ['P'+str(i+1) for i in range(int(num_degen)+1) if i != int(AA_position1)-1 or i != int(AA_position2)-1]
    return [Positions, AA_counts], target_counter_both

def write_output(output_file, AA_target1, AA_position1,AA_target2, AA_position2, common_count, target_count1,target_count2,target_counter_both, num_reads, num_degen, AA_count, Positions):
    '''
    Output information to file
    '''
    if AA_position1 > num_degen:
        output_file.write('Error! Your chosen position is not in the range of target degenerate codons')
    if AA_position2 > num_degen:
        output_file.write('Error! Your chosen position is not in the range of target degenerate codons')
    
    #AA_count_common = AA_count.most_common()
 
    with open(output_file,'w') as output_file:
        output_file.write('The number of times the common, post-degenerate codon sequence appears is '+ str(common_count)+ ' times across '+ str(int(num_reads-1))+ ' reads.\n\n')
        output_file.write('You are interested in '+str(AA_target1)+ ' at position '+str(AA_position1)+'. It appears '+str(target_count1)+ ' times in this file.\n\n')
        output_file.write('You are interested in '+str(AA_target2)+ ' at position '+str(AA_position2)+'. It appears '+str(target_count2)+ ' times in this file.\n\n')
        output_file.write(f"You are interested in {AA_target1} at position {AA_position1} and {AA_target2} at position {AA_position2}. "
                          f"They appear together {target_counter_both} times.\n\n")
        for i in range(len(Positions)-1):
            output_file.write (str(Positions[i]) + '\n') 
            for key, value in AA_count[i].items():
                output_file.write(str(key) + ',' + str(value) + '\n')      

 
        
        
        #for pair in count_full_appearance_sort:
         #   output_file.write (str(pair[0]) + ' ' + str(pair[1]) + '\n')
         
         
if __name__ == '__main__':
    app = QApplication([])
    window = MyWindow()
    app.exec_()
    [NGS_file, common, num_degen, AA_target1, AA_position1, AA_target2, AA_position2, output_file] = user_inputs()
    [degenlist, codons_list, common_count, num_reads] = get_degen_list(NGS_file, common, int(num_degen))
    count_full_AA  = translate_degen_list(degenlist,num_degen) 
    target_count1 = position_counting1 (count_full_AA, AA_target1, AA_position1, num_degen)
    target_count2 = position_counting2 (count_full_AA, AA_target1, AA_position1, num_degen)
    [Positions, AA_counts], target_counter_both = before_after (AA_target1, AA_position1,AA_target2,AA_position2, num_degen, count_full_AA )
    write_output(output_file, AA_target1, AA_position1,AA_target2, AA_position2, common_count, target_count1,target_count2,target_counter_both, num_reads, num_degen, AA_counts, Positions)
    print ("Run completed.")

    

    