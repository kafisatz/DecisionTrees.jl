using PyCall
#check if python has the requird Packages installed

required_python_packages=["pandas","xlsxwriter"] # for testing this functionality you can add some package like "nose" "SymPy" "nltk" "Scapy" ....
try
     for x in required_python_packages
    	PyCall.pywrap(PyCall.pyimport(x))         
     end
 catch
    	@info "One of the required python packages was not found on your system:"
    	println(required_python_packages)
      for x in required_python_packages
        println(x)
      end
    	@info "Trying to install those packages..."
    		# Use eventual proxy info
    		proxy_arg=Array{AbstractString}(undef,0)
    		if haskey(ENV, "http_proxy")
    			push!(proxy_arg, "--proxy")
    			push!(proxy_arg, ENV["http_proxy"])
    		end

    		# Import pip
    		try
    			@pyimport pip
    		catch
    			# If it is not found, install it
    			println("Pip not found on your sytstem. Downloading it.")
    			get_pip = joinpath(dirname(@__FILE__), "get-pip.py")
    			download("https://bootstrap.pypa.io/get-pip.py", get_pip)
    			run(`$(PyCall.python) $(proxy_arg) $get_pip --user`)
    		end

    		println("Installing required python packages using pip")
    		run(`$(PyCall.python) $(proxy_arg) -m pip install --user --upgrade pip setuptools`)
    		run(`$(PyCall.python) $(proxy_arg) -m pip install --user --upgrade pip`)
    		run(`$(PyCall.python) $(proxy_arg) -m pip install --user $(required_python_packages)`)

    	@info "installation of python packages finished."
end
