#Repeat values in structures to match a given size (assuming your structures are mutable)
#This assumes you only want to repeat the array along a single dimension
function RepeatStruct(Structure,Size)

	#Get fields in structure
	FieldsInStruct=fieldnames(typeof(Structure));
	
	for i=1:length(FieldsInStruct)
		#Check field i
		Value=getfield(Structure, FieldsInStruct[i])
		#If field is a float then we repeat it
		if typeof(Value) == Float64 || Int64
			#repeat to predefined size
			Value=repeat([Value],Size[1])
			#Put back in structure
			setfield!(Structure,FieldsInStruct[i],Value)
		else
		#Call the function again for inner structure
			RepeatStruct(Value,Size)
		end
	end

	return Structure

end