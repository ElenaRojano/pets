# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 'gephepred/version'

Gem::Specification.new do |spec|
  spec.name          = "gephepred"
  spec.version       = Gephepred::VERSION
  spec.authors       = ["Elena Rojano, Pedro Seoane"]
  spec.email         = ["elenarojano@uma.es, seoanezonjic@uma.es"]

  spec.summary       = %q{Suite with predictive tools.}
  spec.description   = %q{TODO: Write a longer description or delete this line.}
  spec.homepage      = "TODO: Put your gem's website or public repo URL here."
  spec.license       = "MIT"

  spec.files         = `git ls-files -z`.split("\x0").reject { |f| f.match(%r{^(test|spec|features)/}) }
  spec.bindir        = "bin"
  spec.executables   = spec.files.grep(%r{^bin/}) { |f| File.basename(f) }
  spec.require_paths = ["lib"]

  spec.add_development_dependency "bundler", "~> 1.11"
  spec.add_development_dependency "rake", "~> 10.0"
  spec.add_development_dependency "rspec", "~> 3.0"
  spec.add_dependency "statistics2"
  spec.add_dependency "terminal-table"
  spec.add_dependency "report_html" # ask about this gem

end
