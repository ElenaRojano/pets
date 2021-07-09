# PETS

PETS (Patient Exploration Tools Suite) include three different tools for the analysis of cohorts of patients with pathological phenotypes described in terms of the Human Phenotype Ontology (HPO) and the position their genomic variants clinically determined.

It can (1) determine the quality of information within a patient cohort with Cohort Analyzer (coPatReporter.rb); (2) associate genomic regions with their pathological phenotypes based on the cohort data with Reg2Phen (reg2phen.rb), and (3) predict the possible genetic variants that cause the clinically observed pathological phenotypes using phenotype-genotype association values with Phen2Reg (phen2reg.rb). 

This tool has been developed to be used by the clinical community, to facilitate patient characterisation, help identify where data quality can be improved within a cohort and help diagnose patients with complex disease. Please cite us as Rojano E., Seoane-Zonjic P., Jabato F.M., Perkins J.R., Ranea J.A.G. (2020) Comprehensive Analysis of Patients with Undiagnosed Genetic Diseases Using the Patient Exploration Tools Suite (PETS). In: Rojas I., Valenzuela O., Rojas F., Herrera L., Ortuño F. (eds) Bioinformatics and Biomedical Engineering. IWBBIO 2020. Lecture Notes in Computer Science, vol 12108. Springer, Cham. https://doi.org/10.1007/978-3-030-45385-5_69.


## Installation

Add this line to your application's Gemfile:

```ruby
gem 'pets'
```

And then execute:

    $ bundle

Or install it yourself as:

    $ gem install pets


After installing PETS Gem, R dependencies must be installed. For this, the user must run the following command:
    
    $ install_deps.rb

## Usage

### 1) Cohort Analyzer

Cohort Analyzer measures the phenotyping quality of patient and disease cohorts by calculating multiple statistics to give a general overview of the cohort status in terms of the depth and breadth of phenotyping. It can work with cohorts defined exclusively with HPO terms or with both HPO terms and genomic coordinates.

#### Basic usage of Cohort Analyzer:

We provide an example of use of Cohort Analyzer with a dataset from Vulto-van Silfhout, A.T.; Hehir-Kwa, J.Y.; van Bon, B.W.M.; Schuurs-Hoeijmakers, J.H.M.; Meader, S.; Hellebrekers, C.J.M.; Thoonen, I.J.M.; de Brouwer, A.P.M.; Brunner, H.G.; Webber, C.; Pfundt, R.; de Leeuw, N.; De Vries, B.B.A. Clinical Significance of De Novo and Inherited Copy-Number Variation. Human Mutation 2013, 34, 1679–1687. doi:10.1002/humu.22442.

This dataset includes de novo and inherited CNVs to phenotypes related to intellectual disability/developmental delay occurring alongside multiple congenital anomalies. An example of an input file is available in the example_datasets folder within this repository and the code to execute its analysis is provided below: 

```
coPatReporter.rb -i hummu_congenital_full_dataset.txt -o results -p phenotypes -c chr -d patient_id -s start -e stop -m lin
```

Where: 

- -i -> Input cohort, a tab file with patient identifiers and the list of HPOs characterised for each patient.
- -o -> Output path.
- -p -> Column name with phenotypes.
- -c -> Column name with chromosomes.
- -d -> Column name with patient identifiers.
- -s -> Column name with start genomic coordinate.
- -e -> Column name with final genomic coordinate.
- -m -> Semantic similarity measure method.
- -C -> Maximum number of clusters to display.

Further information with all Cohort Analyzer capabilities for setup can be queried as follows:

```
coPatReporter.rb --help
```

### 2) Reg2Phen

This tool is a search engine that finds phenotypes associated with genomic regions or genes of interest. It uses two input files, one with phenotype-genotype associations previously calculated, and a list of genomic coordinates or gene identifiers to find their HPO associated. We provide an example of use in the example_datasets folder within this repository and the code to execute its analysis is provided below: 

```
reg2phen.rb -t associations_file.txt -p genes.txt -b hpo_file -P -g -H -o results/patient1Genes.txt -F $current/results/patient1Genes.html
```
Where: 

- -t -> Input phenotype-genotype associations file.
- -p -> List of genes to find HPOs associated.
- -b -> HPO obo file.
- -P -> Transform association values in P-values.
- -g -> Set if genes identifiers are provided instead of genome coordinates.
- -H -> Activate HTML reporting.
- -o -> Output folder.
- -F -> Semantic similarity measure method.

Associations between pathological phenotypes and genomic regions provided in this example were calculated with NetAnalyzer (https://rubygems.org/gems/NetAnalyzer, Rojano E. et al (2017). Revealing the Relationship Between Human Genome Regions and Pathological Phenotypes Through Network Analysis. LNCS, 10208:197-207) using randomised DECIPHER data (coordinates in the GRCh37 human genome assembly) and the hypergeometric association method.

### 3) Phen2Reg

Phen2Reg analyses the pathological phenotypes observed in a patient and predicts putative causal genomic regions. As in the case of Reg2Phen, it uses phenotype-genotype associations previously calculated. We provide an example of use in the example_datasets folder within this repository and the code to execute its analysis is provided below:  

```
phen2reg.rb -t associations_file.txt -p example_patient_hpos.txt -i hpo2ci.txt -f hpo_file -T -Q > single_phens.txt
```
Where: 

- -t -> Input phenotype-genotype associations file.
- -p -> List of HPOs characterised for a patient.
- -i -> HPO information coefficients (IC) file.
- -f -> HPO obo file.
- -T -> Deactivate HTML reporting.
- -Q -> Deactivate quality control.

Results are saved in the single_phens.txt output file.

## Development

After checking out the repo, run `bin/setup` to install dependencies. Then, run `rake spec` to run the tests. You can also run `bin/console` for an interactive prompt that will allow you to experiment.

To install this gem onto your local machine, run `bundle exec rake install`. To release a new version, update the version number in `version.rb`, and then run `bundle exec rake release`, which will create a git tag for the version, push git commits and tags, and push the `.gem` file to [rubygems.org](https://rubygems.org).

## Contributing

Bug reports and pull requests are welcome on GitHub at https://bitbucket.org/elenarojano/pets.


## License

The gem is available as open source under the terms of the [MIT License](http://opensource.org/licenses/MIT).

