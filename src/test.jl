struct Person
  myname::String
  age::Int64
end

struct Shoes
   shoesType::String
   colour::String
end

struct Student
   s::Person
   school::String
   shoes::Shoes
end

function printMyActivity(self::Student)
   println("I study at $(self.school) school")
end

struct Employee
   s::Person
   monthlyIncomes::Float64
   company::String
   shoes::Shoes
end

function printMyActivity(self::Employee)
  println("I work at $(self.company) company")
end

gymShoes = Shoes("gym","white")
proShoes = Shoes("classical","brown")

Marc = Student(Person("Marc",15),"Divine School",gymShoes)
MrBrown = Employee(Person("Brown",45),1200.0,"ABC Corporation Inc.", proShoes)

printMyActivity(Marc)
printMyActivity(MrBrown)


Marc = Student(Person("Simone",15),"holy cow",gymShoes)
MrBrown = Employee(Person("White",45),1200.0,"NBC news.", proShoes)

printMyActivity(Marc)
printMyActivity(MrBrown)
