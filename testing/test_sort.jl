using SortingAlgorithms

raw_rel=readdlm("C:\\Users\\bernhard.konig\\Documents\\ASync\\irobot\\data_tmp\\eraw_rel.csv");
orig=deepcopy(raw_rel[:]);
 
r=deepcopy(orig);i=zeros(Int,length(r)); @benchmark mysortperm!(i,r)
r2=deepcopy(orig);i2=zeros(Int,length(r2)); @benchmark sortperm!(i2,r2)
@assert isequal(r2[i2],r[i])

isequal(i,i2)



r=deepcopy(orig); a=sortperm(r)
r=deepcopy(orig); a1,b1=mysortperm(r)

isequal(a,a1)
isequal(raw_rel[a],b1)

r=deepcopy(orig); @btime sort!(r);

r=deepcopy(orig); @btime sort!(r);
r=deepcopy(orig); @btime sort!(r,alg=QuickSort);
r=deepcopy(orig); @btime sort!(r, alg=RadixSort);

r=deepcopy(orig); @btime sortperm(r);
r=deepcopy(orig); @btime mysortperm(r);
r=deepcopy(orig); @btime sortperm(r,alg=QuickSort);
r=deepcopy(orig); @benchmark mysortperm(r)
r=deepcopy(orig); @benchmark sortperm(r,alg=QuickSort)


r=deepcopy(orig);i=zeros(Int,length(r)); @btime sortperm!(i,r,alg=QuickSort);



size(orig)

xx=rand(length(orig));@btime sortperm(xx);
xx=rand(length(orig));@btime sort!(xx);





export mysortperm


function mysortperm!(ii,A)
    const n = length(A)
   #ii = Array(Int64,n)

   for i = 1:n
    ii[i] = i
   end

   B = copy(A)
   quicksort!(B,ii)

   #return ii, B
   return nothing
end # function mysortperm!



function mysortperm(A)
    const n = length(A)
   ii = Array(Int64,n)

   for i = 1:n
    ii[i] = i
   end

   B = copy(A)
   quicksort!(B,ii)

   return ii, B
end # function mysortperm


function quicksort!(A, order, i=1,j=length(A))
# modified from:
# http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#Julia   
    if j > i
          if  j - i <= 50 
              # Insertion sort for small groups is faster than Quicksort
             InsertionSort!(A,order, i,j)
             return A
          end

        #pivot = A[rand(i:j)] # random element of A
        pivot = A[ div(i+j,2) ] 
        left, right = i, j
        while left <= right
            while A[left] < pivot
                left += 1
            end
            while A[right] > pivot
                right -= 1
            end
            if left <= right
                A[left], A[right] = A[right], A[left]
                order[left], order[right] = order[right], order[left]

                left += 1
                right -= 1
            end
        end  # left <= right

        quicksort!(A,order, i,   right)
        quicksort!(A,order, left,j)
    end  # j > i

    return A
end # function quicksort!


function InsertionSort!(A, order, ii=1, jj=length(A))

    for i = ii+1 : jj
        j = i - 1
        temp  = A[i]
        itemp = order[i]

        while true
            if j == ii-1
                 break
            end
            if A[j] <= temp
                 break
            end
            A[j+1] = A[j]
            order[j+1] = order[j]
            j -= 1
        end

        A[j+1] = temp
        order[j+1] = itemp
    end  # i

return
end # function InsertionSort!