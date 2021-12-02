require 'numo/narray'
require 'npy'

# TODO: Ver de llamar a Netanalyzer aqui para usar su expansion del Numo Narray para meter metodos para tranformar un hash en matriz y viceversa. 
# Tb meterle metodos para escribir y lerr matrices con npy
def hash2bin_matrix(hash)
  x_names_indx = get_hash_values_idx(hash)
  y_names = hash.keys
  x_names = x_names_indx.keys
   # row (y), cols (x)
  matrix = Numo::DFloat.zeros(hash.length, x_names.length)
  i = 0
  hash.each do |id, items|
    items.each do |item_id|
      matrix[i, x_names_indx[item_id]] = 1
    end
    i += 1
  end
  return matrix, y_names, x_names
end

def get_hash_values_idx(hash)
  x_names_indx = {}
  i = 0
  hash.each do |k, items|
    items.each do |item_id|
      query = x_names_indx[item_id]
      if query.nil?
        x_names_indx[item_id] = i
        i += 1
      end
    end
  end
  return x_names_indx
end

def hash2weighted_matrix(profiles_similarity)
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

def save_matrix(matrix, matrix_filename, x_axis_names, x_axis_file, y_axis_names=nil, y_axis_file=nil)
  File.open(x_axis_file, 'w'){|f| f.print x_axis_names.join("\n") }
  File.open(y_axis_file, 'w'){|f| f.print y_axis_names.join("\n") } if !y_axis_names.nil?
  Npy.save(matrix_filename, matrix)
end
