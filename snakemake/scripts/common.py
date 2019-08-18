
import logging
from logging.handlers import FileHandler

class BufferedOutputHandler:
	buffer_size = 4000

	def __init__(self, output_file, log_file, output_type='fasta'):
		self.output_file = output_file
		self.output_handle = open(output_file, 'a')
		self.output_type = output_type
		self.buffer = []
		self.total_count = 0
		self.my_logger = LoggerUtility('BufferedOutputHandler_%s' % output_file, log_file).get_logger()

		self.my_logger.debug('Opened up the buffered file')

	def add_record(self, record):
		"""
		Adds a record to the buffer. If the number of records exceeds the buffer size,
		the records are written out to save memory.
		:param record: The record to be added
		:return: None
		"""
		self.buffer.append(record)
		if len(self.buffer) >= self.buffer_size:
			count = SeqIO.write(self.buffer, self.output_handle, self.output_type)
			self.output_handle.flush()
			self.total_count += self.buffer_size
			self.buffer = []
			self.my_logger.debug('Buffering %i reads (total %i)' % (count, self.total_count))

	def close_out(self, zip_output=True):
		"""
		Closes out the handle to write out records
		:param zip_output: Whether or not to zip the output
		:return: None
		"""
		SeqIO.write(self.buffer, self.output_handle, self.output_type)
		self.output_handle.flush()
		self.output_handle.close()
		self.my_logger.debug('Closing out (%i reads added, total %i)' % (len(self.buffer),
		                                                                 self.total_count + len(self.buffer)))

class LoggerUtility:

	def __init__(self, logger_name, log_file, level=logging.DEBUG):
		self.FORMATTER = logging.Formatter("%(asctime)s — %(name)s — %(levelname)s — %(message)s")
		self.logger = logging.getLogger(logger_name)
		self.logger.setLevel(level)  # better to have too much log than not enough
		self.logger.addHandler(self.get_console_handler())
		self.logger.addHandler(self.get_file_handler(log_file))
		self.logger.propagate = False

	def get_console_handler(self):
		"""
		Creates a logger to go to the console
		:return: console handle
		"""
		console_handler = logging.StreamHandler(sys.stdout)
		console_handler.setFormatter(self.FORMATTER)
		return console_handler

	def get_file_handler(self, log_file):
		"""
		Creates a file handler to write to a file
		:param log_file: The file to write to
		:return: file handler
		"""
		file_handler = FileHandler(log_file)
		file_handler.setFormatter(self.FORMATTER)
		return file_handler

	def get_logger(self):
		return self.logger
