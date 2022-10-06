# Needs define ROOT_PATH constant in file requiring this file
COMMON_OPTPARSE = File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'pets', 'common_optparse.rb'))
REPORT_FOLDER = File.expand_path(File.join(ROOT_PATH, '..', 'templates'))
EXTERNAL_DATA = File.expand_path(File.join(ROOT_PATH, '..', 'external_data'))
EXTERNAL_CODE = File.expand_path(File.join(ROOT_PATH, '..', 'external_code'))
HPO_FILE = File.join(EXTERNAL_DATA, 'hp.json')
MONDO_FILE = File.join(EXTERNAL_DATA, 'mondo.obo')
IC_FILE = File.join(EXTERNAL_DATA, 'uniq_hpo_with_CI.txt')
