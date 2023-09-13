using LinearAlgebra
using DelimitedFiles

include("helpers.jl")


#main functions 



function alg1main()
    
    println("Input a matrix of integers. Press Ctrl D (LINUX) or Ctrl Z (WINDOWS) when done.")
    m1 = readdlm(stdin)
    if validShape(m1) != -1
        m1 = Int.(m1)

        if valid(m1)
            n = size(m1,1)
            zarreMatrix = ZARREform(m1, n) 
            upperT = upperTriangularize(zarreMatrix, n) 
            reconstructablePrinter(reconstructable(upperT,n))   
        else
            errorMessage(m1, n) 
        end 
    end   

end


function alg2main()
    print("How many induced orders will you input: ")
    n = parse(Int,readline())
    println(" ")

    induced = ZARREList(inducedMatrices(Int.(matrix9),n),n)

    println("INDUCED: ",induced)
    reconstructedMatrix = reconstruction(induced,n)
    println(reconstructedMatrix)


end

function valid(m1) #determines whether an nxn matrix is valid 
    if positiveFirstCol(m1) && (round(det(m1), RoundToZero)  != 0.0) 
        return true
    else
        return false
    end
end

function validShape(m1) #determines whether a matrix is of nxn form. Returns n, or -1 if invalid. 
    if typeof(m1) != Matrix{Float64} || size(m1,1) != size(m1,2)
        println("INVALID INPUT - FAILED SHAPE VALIDATION")
        return -1
    else
        return size(m1,1)
    end 
end


function ZARREform(m,n) #n represents the size of the matrix.
    for i in 1:n-1
        c = nonZeroIdx(m[i,:])
        a = nonZeroVal(m[i,:])
        for j in i+1:n
            if m[j,c] != 0
                m[j,:] = addOrSubtractRows(rowMultiplication(m[j,:],abs(a)),rowMultiplication(m[i,:],(a * m[j,c])/abs(a)))
                d = gcd(m[j,:])
                println(d)
                m[j,:] = rowMultiplication(m[j,:], 1/d)
            end
        end
    end
    return m
end



function upperTriangularize(zarre, n) #creates the upper triangularized matrix by permuting the appropriate columns (go over)
    r = 0
    newMatrix = zeros(n,n)
    for origRow in eachrow(zarre)
        r += 1
        newCol = nonZeroIdx(origRow) #newCol is index of leading unit in origRow 
        newMatrix[:,r] = eachcol(zarre)[newCol]
    end
    return newMatrix
end


function reconstructable(A,n)
    reconstructable = true

    for i in 1:n-2 #for each row, starting with one and ending with n-2    
        B = A[:, 1:end .!= i] #slice out the i'th column and save the rest to a new matrix 
        if nonZeroIdx(B[i,:]) == 0 
            return reconstructable
        else
            if n > 3
                k = nonZeroIdx(B[i,:]) #save the first non-zero column in row i
                rowKPlusOne = B[k+1, :]

                c = B[k+1,k] / B[i,k]
                println(c)

                for t in k+1:n #k < t <= n
                    if B[k+1][t] == c * B[i,t] #if multiple
                            continue
                    else
                            reconstructable = false
                    end
                end
            else
                reconstructable = false # n=3 and first row isn't all zeros 
            end

        end
    end
    return reconstructable 
end


function input(n) #shape of reconstructed matrix
    println("Enter ", n, " ", n-1,"x",n-1, " induced order matrices— after each row, press enter. ")
    v = Matrix{Int64}[]
    for i in 1:n
        println("Begin induced order ", i, ": ")
        m1 = zeros(n-1,n-1)
        for j in 1:n-1
            m1[j,:] = parser(readline())
        end
        push!(v,Int.(m1))
        println(" ")
        println("You have entered: ", m1,)
        println(" ")

        if valid(m1)
            continue
        else
            return []
        end
    end
    println(" ")
    println("You have entered the following induced orders:")
    
    ct = 0
    for matrix in v
        ct += 1
        println(ct, " induced order: ", matrix)
    end
    return v
end


function reconstruction(inducedOrdersList, n) #n is the size of the recsonstructed matrix 
    m = zeros(Int, n,n)
    secondToLastInduced = inducedOrdersList[n-1]
    lastInduced = inducedOrdersList[n]

    secondToLastwithoutLastCol = secondToLastInduced[:, 1:end .!= n-1]
    lastwithoutLastCol = lastInduced[:, 1:end .!= n-1]

    if rebuildable(secondToLastwithoutLastCol, lastwithoutLastCol)
        println("THERE IS AT LEAST 1 MATRIX WITH THESE INDUCED ORDERS...")
        rozsl = zeroRowInMatrix(secondToLastwithoutLastCol)
        rozl = zeroRowInMatrix(lastwithoutLastCol)

        if rozsl < rozl
            println("case 1")
            m[rozsl, n] = secondToLastInduced[rozsl, n-1] #rest of entries in row are already zero 
            m[rozl+1,n-1] = lastInduced[rozl, n-1]

            for j in 1:n
                row = m[j,:]
              

                if j < rozsl
                    fnzsl = nonZeroVal(secondToLastInduced[j,:])
                    fnzl = nonZeroVal(lastInduced[j,:])

                    for k in 1:n-1
                        m[j,k] = (abs(fnzsl)/gcd(abs(fnzsl),abs(fnzl))) * lastInduced[j,k]
                    end
                    m[j,n] = abs(fnzl)/gcd(abs(fnzsl),abs(fnzl)) * secondToLastInduced[j,n-1]
                end

                if j > rozsl && j < rozl+1
                    fnzsl = nonZeroVal(secondToLastInduced[j,:])
                    fnzl = nonZeroVal(lastInduced[j-1,:])
                    for k in 1:n-1
                        m[j,k] = (abs(fnzsl)/gcd(abs(fnzsl),abs(fnzl))) * lastInduced[j-1,k]
                    end
                    m[j,n] = abs(fnzl)/gcd(abs(fnzsl),abs(fnzl)) * secondToLastInduced[j,n-1]
                end


                if j > rozl+1 
                    fnzsl = nonZeroVal(secondToLastInduced[j-1,:])
                    fnzl = nonZeroVal(lastInduced[j-1,:])
                    for k in 1:n-1
                        m[j,k] = (abs(fnzsl)/gcd(abs(fnzsl),abs(fnzl))) * lastInduced[j-1,k]
                    end
                    m[j,n] = abs(fnzl)/gcd(abs(fnzsl),abs(fnzl)) * secondToLastInduced[j-1,n-1]         
                end

            end


  

        elseif rozsl > rozl
            println("case 2")
            m[rozsl+1, n] = secondToLastInduced[rozsl, n-1] #flipped here 'rozsl+1' instead of 'rozl+1'
            m[rozl,n-1] = lastInduced[rozl, n-1]
   

            for j in 1:n
                
                if j < rozl
                    fnzsl = nonZeroVal(secondToLastInduced[j,:])
                    fnzl = nonZeroVal(lastInduced[j,:])

                    for k in 1:n-1
                        m[j,k] = (abs(fnzsl)/gcd(abs(fnzsl),abs(fnzl))) * lastInduced[j,k]
                    end
                    m[j,n] = abs(fnzl)/gcd(abs(fnzsl),abs(fnzl)) * secondToLastInduced[j,n-1]


                elseif j > rozl && j < rozsl + 1
                    fnzsl = nonZeroVal(secondToLastInduced[j-1,:]) #flipped here 
                    fnzl = nonZeroVal(lastInduced[j,:])

                    for k in 1:n-1
                        m[j,k] = (abs(fnzsl)/gcd(abs(fnzsl),abs(fnzl))) * lastInduced[j,k]
                    end
                    m[j,n] = abs(fnzl)/gcd(abs(fnzsl),abs(fnzl)) * secondToLastInduced[j-1,n-1]

                elseif j > rozsl + 1
                    fnzsl = nonZeroVal(secondToLastInduced[j-1,:])
                    fnzl = nonZeroVal(lastInduced[j-1,:])

                    for k in 1:n-1
                        m[j,k] = (abs(fnzsl)/gcd(abs(fnzsl),abs(fnzl))) * lastInduced[j-1,k]
                    end
                    m[j,n] = abs(fnzl)/gcd(abs(fnzsl),abs(fnzl)) * secondToLastInduced[j-1,n-1]


                end
            end


        elseif rozsl == rozl && rozl != n-1
            println("case 3")

            s = nonZeroIdx(lastInduced[rozl+1, :])
            sInduced = inducedOrdersList[s]
            m[rozl, n-1] = sInduced[rozl,n-2]
            m[rozl, n] = sInduced[rozl, n-1]

            println("m right before loop: ", m)

            for j in 1:rozl-1
                fnzl = nonZeroVal(lastInduced[j,:])
                fnzsl = nonZeroVal(secondToLastInduced[j,:])

                for k in 1:n-1
                    m[j,k] = abs(fnzsl)/gcd(abs(fnzl),abs(fnzsl)) * lastInduced[j,k]
                end
                m[j,n] = abs(fnzl)/gcd(abs(fnzl),abs(fnzsl)) * secondToLastInduced[j,n-1]
            end

            k = rozl + 1
            match = false  

            while (k <= n && match == false)

                if sInduced[rozl,n-2] != 0
                    m[k,n] = 1             
                elseif sInduced[rozl, n-2] == 0
                    m[k,n-1] = 1       
                end 

                for j in rozl+1:k-1
                    fnzl = nonZeroVal(lastInduced[j-1,:])
                    fnzsl = nonZeroVal(secondToLastInduced[j,:])

                    for p in 1:n-1
                        m[j,p] = abs(fnzsl)/gcd(abs(fnzl),abs(fnzsl)) * lastInduced[j-1,p]
                    end
                    m[j,n] = abs(fnzl)/gcd(abs(fnzl),abs(fnzsl)) * secondToLastInduced[j,n-1]
                end

                for j in k+1:n
                    fnzl = nonZeroVal(lastInduced[j-1,:])
                    fnzsl = nonZeroVal(secondToLastInduced[j-1,:])

                    for p in 1:n-1
                        m[j,p] = abs(fnzsl)/gcd(abs(fnzl),abs(fnzsl)) * lastInduced[j-1,p]
                    end
                    m[j,n] = abs(fnzl)/gcd(abs(fnzl),abs(fnzsl)) * secondToLastInduced[j-1,n-1]


                end
                println("finished the case")
                println("m after case: ", m)
                inducedM = inducedMatrices(m,n)
                
                if inducedM == inducedOrdersList
                    match = true
                
                else
                    m[k:end,:] = zeros(Int,n-k,n)
                    k += 1
                end
            end  
                
            
        elseif rozsl == rozl && rozl == n-1
            println("case 4 - unable to reconstruct matrix.")
        end




    else
        println("THESE INDUCED ORDERS DO NOT CORRESPOND TO A VALID MATRIX")
    end

    return m

end


function inducedMatrices(matrix, n) 
    inducedMatricesList = Matrix{Int64}[]
    if valid(matrix) #only continue if the matrix is valid (valid() function is from algorithm 1) COME BACK

    C = cMaker(matrix, n) #create the C matrix 
    for i in 1:n
        A = matrix[:, 1:end .!= i ] #eliminate the i'th column

        s = C[2, i]
        a = C[4, i]

        if a!= 0
            e = C[2, C[5, i]]
            b = C[3,C[5,i]]
        end


        while  a!= 0  
            A[e, :] = addOrSubtractRows(rowMultiplication(A[e,:], abs(a)),rowMultiplication(rowMultiplication(A[s,:], abs(b)), (a*b)/abs(a*b)))
            s = e
            a = nonZeroVal(A[s,:]) #returns the value
            c = nonZeroIdx(A[s,:])
            e = C[2, c+1]
            b = C[3, c+1]     
        end

        final = A[1:end .!= s, :]
        
        for i in 1:n-1
            d = gcd(final[i,:])
            final[i,:] = rowMultiplication(final[i,:], 1/d )
        end

        push!(inducedMatricesList, final)
    end

    else

        errorMessage(matrix, n) #if it's not valid, print out an error message (errorMessage() function is from algorithm 1)

    end
    return inducedMatricesList

end












