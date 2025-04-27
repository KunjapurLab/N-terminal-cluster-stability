from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QLabel, QLineEdit, QPushButton, QFileDialog
from Bio import SeqIO

class MyWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle('FASTQ to FASTA Converter')
        layout = QVBoxLayout()

        # File browser section for FASTQ file
        self.file_browser_label = QLabel('FASTQ Data File:')
        self.file_browser_input = QLineEdit()
        self.file_browser_button = QPushButton('Browse for FASTQ file')
        self.file_browser_button.clicked.connect(self.browse_file)

        layout.addWidget(self.file_browser_label)
        layout.addWidget(self.file_browser_input)
        layout.addWidget(self.file_browser_button)

        # Output file name input
        self.output_label = QLabel('Output FASTA File Name:')
        self.output_input = QLineEdit()
        layout.addWidget(self.output_label)
        layout.addWidget(self.output_input)

        # Button to submit and convert
        self.button_submit = QPushButton('Convert to FASTA')
        self.button_submit.clicked.connect(self.convert_fastq_to_fasta)
        layout.addWidget(self.button_submit)

        self.setLayout(layout)
        self.show()

    def browse_file(self):
        file_dialog = QFileDialog()
        if file_dialog.exec_():
            file_path = file_dialog.selectedFiles()[0]
            self.file_browser_input.setText(file_path)

    def convert_fastq_to_fasta(self):
        fastq_file = self.file_browser_input.text()
        output_file = self.output_input.text()

        if fastq_file and output_file:
            output_file = output_file if output_file.endswith('.fasta') else output_file + '.fasta'
            self.fastq_to_fasta(fastq_file, output_file)
            print(f"Conversion completed. Output saved to {output_file}")
        else:
            print("Please provide both input and output file names.")

    def fastq_to_fasta(self, input_file, output_file):
        # Open the FASTQ file directly without gzip, as it's uncompressed
        with open(input_file, 'r') as fastq:
            with open(output_file, 'w') as fasta_out:
                for record in SeqIO.parse(fastq, 'fastq'):
                    seq = str(record.seq).replace('\n', '')
                    fasta_out.write(f'>{record.id}\n{seq}\n')

if __name__ == '__main__':
    app = QApplication([])
    window = MyWindow()
    app.exec_()
