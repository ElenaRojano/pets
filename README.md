# Pets

Pets (Patient exploration tools suite) include tools for the analysis of cohorts of patients with pathological phenotypes described in terms of the Human Phenotype Ontology (HPO) and the position their genomic variants clinically determined. 

Pets include tools to (1) perform cohort analysis (coPatReporter.rb), (2) searching for pathological phenotypes associated with a genomic region of interest (reg2phen.rb) and (3) predict regions of the genome that potentially lead to the pathological phenotypes observed in a patient (phen2reg.rb).

Associations between pathological phenotypes and genomic regions (using genomic coordinates from GRCh37 human assembly) are previously calculated using NetAnalyzer (https://rubygems.org/gems/NetAnalyzer). Please cite us as Rojano E. et al (2017). Revealing the Relationship Between Human Genome Regions and Pathological Phenotypes Through Network Analysis. LNCS, 10208:197-207.

## Installation

Add this line to your application's Gemfile:

```ruby
gem 'pets'
```

And then execute:

    $ bundle

Or install it yourself as:

    $ gem install pets


After installing PETS Gem, R dependencies must be installed by running the script:

    $ path_to_gem/external_code/install_R_dependencies.R

## Usage

TODO: Write usage instructions here

## Development

After checking out the repo, run `bin/setup` to install dependencies. Then, run `rake spec` to run the tests. You can also run `bin/console` for an interactive prompt that will allow you to experiment.

To install this gem onto your local machine, run `bundle exec rake install`. To release a new version, update the version number in `version.rb`, and then run `bundle exec rake release`, which will create a git tag for the version, push git commits and tags, and push the `.gem` file to [rubygems.org](https://rubygems.org).

## Contributing

Bug reports and pull requests are welcome on GitHub at https://bitbucket.org/elenarojano/pets.


## License

The gem is available as open source under the terms of the [MIT License](http://opensource.org/licenses/MIT).

