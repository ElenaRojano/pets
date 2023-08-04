#! /usr/bin/env ruby

ROOT_PATH = File.dirname(__FILE__)
EXTERNAL_CODE = File.expand_path(File.join(ROOT_PATH, '..', 'external_code'))
$: << File.expand_path(File.join(ROOT_PATH, '..', 'lib', 'pets'))
require 'pets'

system_call(EXTERNAL_CODE, 'install_R_dependencies.R', '')