$:.unshift File.join(File.dirname(__FILE__), ".")
require 'libsbmlsim'

res = Libsbmlsim::simulateSBMLFromFile('../../hogeMAPK.xml', 4000.0, 0.1, 100, 1, Libsbmlsim::MTHD_RUNGE_KUTTA, 0)
puts res.getErrorMessageAtIndex(0) if res.getNumOfErrors > 0
res = Libsbmlsim::simulateSBMLFromFile('../../MAPK.xml', 4000.0, 0.1, 100, 1, Libsbmlsim::MTHD_RUNGE_KUTTA, 0)

numOfRows = res.getNumOfRows
puts "numOfRows: #{numOfRows}"

numOfSp = res.getNumOfSpecies
puts "numOfSpeies: #{numOfSp}"

numOfParam = res.getNumOfParameters
puts "numOfParameters: #{numOfParam}"

numOfComp = res.getNumOfCompartments
puts "numOfCompartments: #{numOfComp}"

timeName = res.getTimeName
puts "timeName: #{timeName}"

puts "Species Name"
numOfSp.times do |n|
  sname = res.getSpeciesNameAtIndex(n)
  puts "  #{sname}"
end

puts "Parameter Name"
numOfParam.times do |n|
  pname = res.getParameterNameAtIndex(n)
  puts "  #{pname}"
end

puts "Compartment Name"
numOfComp.times do |n|
  cname = res.getCompartmentNameAtIndex(n)
  puts "  #{cname}"
end

puts "Values"
numOfRows.times do |i|
  t = res.getTimeValueAtIndex(i)
  print "#{t}"
  numOfSp.times do |j|
    sname = res.getSpeciesNameAtIndex(j)
    v = res.getSpeciesValueAtIndex(sname, i)
    print " #{v}"
  end
  print "\n"
  break if i == 10;
end

