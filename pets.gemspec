
lib = File.expand_path("../lib", __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require "pets/version"

Gem::Specification.new do |spec|
  spec.name          = "pets"
  spec.version       = Pets::VERSION
  spec.authors       = ["Elena Rojano, Pedro Seoane"]
  spec.email         = ["elenarojano@uma.es, seoanezonjic@uma.es"]

  spec.summary       = %q{Suite with predictive tools.}
  spec.description   = %q{PETS suite includes three different tools. CohortAnalyzer performs the calculation of several statistics that gives an overview of a cohort of patients to analyse. Reg2Phen uses associations between pathological phenotypes and regions of the genome (these associations can be calculated from the cohort of patients if they include genotypic & phenotypic information using NetAnalyzer, another Ruby gem) to find, for a given genomic region, which pathological phenotypes have been associated with that region. The third tool, Phen2Reg, is a predictor that using the same associations as Reg2Phen, predicts which genomic regions can be the cause of a list of pathological phenotypes observed in a patient.}
  spec.homepage      = "https://bitbucket.org/elenarojano/reg2phen/src/master/bin/reg2phen.rb"
  spec.license       = "MIT"

  # Prevent pushing this gem to RubyGems.org. To allow pushes either set the 'allowed_push_host'
  # to allow pushing to a single host or delete this section to allow pushing to any host.
#  if spec.respond_to?(:metadata)
#    spec.metadata["allowed_push_host"] = "TODO: Set to 'http://mygemserver.com'"
#
#    spec.metadata["homepage_uri"] = spec.homepage
#    spec.metadata["source_code_uri"] = "TODO: Put your gem's public repo URL here."
#    spec.metadata["changelog_uri"] = "TODO: Put your gem's CHANGELOG.md URL here."
#  else
#    raise "RubyGems 2.0 or newer is required to protect against " \
#      "public gem pushes."
#  end

  # Specify which files should be added to the gem when it is released.
  # The `git ls-files -z` loads the files in the RubyGem that have been added into git.
  spec.files         = Dir.chdir(File.expand_path('..', __FILE__)) do
    `git ls-files -z`.split("\x0").reject { |f| f.match(%r{^(test|spec|features)/}) }
  end
  spec.bindir        = "bin"
  spec.executables   = spec.files.grep(%r{^bin/}) { |f| File.basename(f) }
  spec.require_paths = ["lib"]

  spec.add_development_dependency "bundler", "~> 2.0"
  spec.add_development_dependency "rake", "~> 10.0"
  spec.add_development_dependency "rspec", "~> 3.0"
  spec.add_dependency "statistics2"
  spec.add_dependency "terminal-table"
  spec.add_dependency "semtools", "~> 0.1.0" 
  spec.add_dependency "report_html" 
  spec.add_dependency "numo-narray" 
  spec.add_dependency "npy" 
  spec.add_dependency "parallel"  

end

