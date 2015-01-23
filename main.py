import sys
import os
import os.path
import datetime
import shutil
import utils
import logging
import argparse
from utils import STEPS

class State(object):
  base_dir = os.path.dirname(os.path.realpath(__file__))
  working_dir = None
  
  def set_working_dir(self, directory):
    self.working_dir = directory
    self.log = os.path.join(self.working_dir, 'log.log')
    self.last_state = os.path.join(self.working_dir, 'state')
    self.input_pdb = os.path.join(self.working_dir, 'input.pdb')
    self.reference_sequence = os.path.join(self.working_dir, 'reference.fasta')
    self.distances_matrix = os.path.join(self.working_dir, 'distances_matrix.txt')
    self.blast1xml = os.path.join(self.working_dir, 'blast1.xml')
    self.blast1fasta = os.path.join(self.working_dir, 'blast1.fasta')
    self.blast1aln = os.path.join(self.working_dir, 'blast1.aln')
    self.blast2xml = os.path.join(self.working_dir, 'blast2.xml')
    self.blast2fasta = os.path.join(self.working_dir, 'blast2.fasta')
    self.blast2aln = os.path.join(self.working_dir, 'blast2.aln')
    self.blast3xml = os.path.join(self.working_dir, 'blast3.xml')
    self.blast3fasta = os.path.join(self.working_dir, 'blast3.fasta')    
    self.msa = os.path.join(self.working_dir, 'msa.fasta')
    self.correlations = os.path.join(self.working_dir, 'correlations.txt')
    self.correlations_reduced = os.path.join(self.working_dir, 'correlations_reduced.txt')
    self.correlation_for_pair = os.path.join(self.working_dir, 'correlation_for_pair.txt')
    
  def __getitem__(self, name):
    return self.__getattribute__(name)
    
  
class CmdLineParser(object):
  parsed = None
  def __init__(self):
    self.parser = argparse.ArgumentParser(description='bioinformatics project tool')    
    subparsers = self.parser.add_subparsers(help='sub-command help')
    
    parser_new = subparsers.add_parser('new', help='start process a new sequence')
    parser_new.set_defaults(mode='new')
    parser_new.add_argument('-t', '--to_step', type=int, default='999', help="run till this step")
    parser_new.add_argument('input', type=str, help="input pdb file")  
    
    parser_continue = subparsers.add_parser('continue', help='continue from previous run')
    parser_continue.set_defaults(mode='continue')
    parser_continue.add_argument('-f', '--from_step', type=int, help="starting at this step. may be constant (1, 2, 3, ...) or 0 to continue after last successull step")  
    parser_continue.add_argument('-t', '--to_step', type=int, default='999', help="run till this step")    
    parser_continue.add_argument('results_directory', nargs='?', help="results directory")
  
  def parse(self, state):
    if not self.parsed:
      self.parsed = self.parser.parse_args()
    if self.parsed.mode == 'new':
      self.mode = 'new'
      self.input = self.parsed.input
      self.from_step = 1
      self.to_step = min(int(self.parsed.to_step), len(STEPS))
    elif self.parsed.mode == 'continue':
      self.mode = 'continue'
      if self.parsed.results_directory:      
        self.results_dir = os.path.abspath(self.parsed.results_directory)
      else:
        # find latest results
        master_results_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'results')
        results = os.listdir(master_results_dir)
        results.sort()
        self.results_dir = os.path.join(master_results_dir, results[-1])
        logging.info('continue on latest results directory - %s' % (self.results_dir))
      state.set_working_dir(self.results_dir)
      if self.parsed.from_step:
        self.from_step = self.parsed.from_step
      else:
        self.from_step = get_last_successfull_step(state)+1
        logging.info('continue from after successull step - %d' % (self.from_step))
      self.to_step = min(int(self.parsed.to_step), len(STEPS))

def get_last_successfull_step(state):
  file_path = state.last_state
  with file(file_path, 'r') as f:
    return int(f.read())
    
def set_last_successfull_step(state, step):
  file_path = state.last_state
  with file(file_path, 'w') as f:
    return f.write(str(step))
  
  
def import_new_sequence(state, input_file_path):
  # create output directory
  output_folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'results', str(datetime.datetime.now().strftime("%Y-%m-%d %H-%M-%S")))
  os.mkdir(output_folder)
  state.set_working_dir(output_folder)
  set_last_successfull_step(state, 0)
  # copy input file
  shutil.copyfile(input_file_path, state.input_pdb)

def set_console_logger():
  # set logger
  logging.basicConfig(level=logging.DEBUG, format='%(asctime)-12s: %(levelname)-8s %(message)s', datefmt='%H:%M:%S')
  logging.debug('program started')

def set_file_logger(state):
  handler = logging.FileHandler(state.log, mode='a')
  handler.setLevel(logging.DEBUG)
  formatter = logging.Formatter('%(asctime)-12s: %(levelname)-8s %(message)s',  datefmt='%d/%m/%Y %H:%M:%S')
  handler.setFormatter(formatter)
  logging.getLogger('').addHandler(handler)
  logging.debug('file logger attached')

def main(argv):
  set_console_logger()
  
  # Init
  state = State()
  
  # Parse
  parser = CmdLineParser()
  parser.parse(state)
  
  # Preparations
  if parser.mode == 'new':
    # create a new folder structure
    import_new_sequence(state, parser.input)

  # Attach to file logger
  set_file_logger(state)

  if parser.from_step > parser.to_step:
    logging.info('no steps to run')
    return
    
  # Clear prev step results
  if parser.mode == 'continue':
    logging.info('deleting previous results (steps %d-%d)' % (parser.from_step, len(STEPS)))
    set_last_successfull_step(state, parser.from_step-1)
    for step in xrange(parser.from_step, len(STEPS)+1):
      current_step_type = STEPS[step-1]    
      current_step = current_step_type()      
      current_step.clean(state)
      
  # Run Steps
  for step in xrange(parser.from_step, parser.to_step+1):  
    current_step_type = STEPS[step-1]    
    current_step = current_step_type()
    logging.info('running step %d - %s' % (step, current_step.name))
    
    current_step.run(state)
    set_last_successfull_step(state, step)

if __name__=="__main__":
  main(sys.argv)