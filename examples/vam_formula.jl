using VirtualAgeModels
using DataFrames

df2 = DataFrame(x=1, y=[1,2,3, 4])



macro assign(name, e)
    Expr(:(=), esc(name), e)
end

@assign a df2

a

ex = Expr(:(=), :a, 1)
eval(ex)
a

ex2 = Meta.parse("a=1")
ex2.head

dump(ex2)

Meta.show_sexpr(ex2)

macro assign2(name, e)
        :($esc(name) = e)
end

df2

@assign2 a2 df2
a2