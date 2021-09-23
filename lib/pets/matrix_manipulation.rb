# TODO: Ver de llamar a Netanalyzer aqui para usar su expansion del Numo Narray para meter metodos para tranformar un hash en matriz y viceversa. 
# Tb meterle metodos para escribir y lerr matrices con npy

def generate_patient_hpo_matrix(patient_data, cohort_hpos)
  matrix = []
  n = cohort_hpos.length
  patient_data.each do |pat_id, pat_hpos|
    vector = Array.new(n, 0)
    pat_hpos.each do |hpo|
      vector[cohort_hpos.index(hpo)] = 1
    end
    matrix << vector
  end
  return matrix
end

def generate_patient_hpo_matrix_numo(patient_data, cohort_hpos)
  y_names = patient_data.keys
  x_names = cohort_hpos
  x_names_indx = {}
  cohort_hpos.each_with_index{|hp,x| x_names_indx[hp]=x}
  # row (y), cols (x)
  matrix = Numo::DFloat.zeros(patient_data.length, cohort_hpos.length)
  i = 0
  patient_data.each do |pat_id, pat_hpos|
    pat_hpos.each do |hp|
      matrix[i, x_names_indx[hp]] = 1
    end
    i += 1
  end
  return matrix, y_names, x_names
end

def format_profiles_similarity_data(profiles_similarity)
  matrix = []
  element_names = profiles_similarity.keys
  matrix << element_names
  profiles_similarity.each do |elementA, relations|
    row = [elementA]
    element_names.each do |elementB|
      if elementA == elementB
        row << 'NA'
      else
        query = relations[elementB]
        if !query.nil?
          row << query
        else
          row << profiles_similarity[elementB][elementA]
        end
      end
    end
    matrix << row
  end
  matrix[0].unshift('pat')
  return matrix
end


def format_profiles_similarity_data_numo(profiles_similarity)
  element_names = profiles_similarity.keys
  matrix = Numo::DFloat.zeros(element_names.length, element_names.length)
  i = 0
  profiles_similarity.each do |elementA, relations|
    element_names.each_with_index do |elementB, j|
      if elementA != elementB
        query = relations[elementB]
        if !query.nil?
          matrix[i, j] = query
        else
          matrix[i, j] = profiles_similarity[elementB][elementA]
        end
      end
    end
    i += 1
  end
  return matrix, element_names
end

def read_excluded_hpo_file(file)
  excluded_hpo = []
  File.open(file).each do |line|
    excluded_hpo << line.chomp
  end
  return excluded_hpo
end