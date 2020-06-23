function myfunc()
    A = rand(200, 200, 400)
    maximum(A)
end

myfunc()

using Profile

@profile myfunc()

using ProfileView

Profile.clear()  # in case we have any previous profiling data
@profile myfunc()
ProfileView.view()