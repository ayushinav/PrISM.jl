using Pkg
Pkg.activate(".")

using Distributions
import Distributions:Distribution
using MT
import MT:boltz_k, gas_R, arrh_dry, arrh_wet


abstract type AbstractRockphyModel end
abstract type AbstractRockphyModelDistribution end


mutable struct Dai_Karato2014{F1, F2} <: AbstractCondModel
    T::F1
    Ch2o_ol::F2
end

params_Dai_Karato2014 = (;
    r = 1f0,
    Cw = 460f0,
    s100 = 10.0f0^74.0f-2,
    s010 = 10.0f0^51.0f-2,
    s001 = 10.0f0^30.0f-2,
    h100 = 75.0f0,
    h010 = 73.0f0,
    h001 = 71.0f0,
    sw100 = 10.0f0^4.52f0,
    sw010 = 10.0f0^1.99f0,
    sw001 = 10.0f0^1.08f0,
    hw100 = 140.0f0,
    hw010 = 101.0f0,
    hw001 = 87.0f0
)

default_params_SEO3 = deepcopy(params_SEO3)

function forward(m::Dai_Karato2014, p, params = default_params_Dai_Karato2014)
    @unpack s = params

    σ_001 = 
    s100.*exp(-h100./(R.*T)) + sw100.*exp(-hw100./(R.*T))


end


mutable struct Zhang2012{F1, F2} <: AbstractCondModel
    T::F1
    Ch2o_opx::F2
end

params_Zhang2012 = (
    S_pol = 3.99f0,
    H_pol = 1.88f0,
    S_hyd = 10f0 ^2.58f0,
    H_hyd = 84f-2,
    a = 8f-2,
)

default_params_Zhang_2012 = deepcopy(params_Zhang2012)

function forward(m::Zhang2012, p, params = default_params_Zhang2012)
    @unpack S_pol, H_pol, S_hyd, H_hyd, a = params
    
    σ_pol =  @. arrh_dry(S_pol, H_pol, boltz_k, m.T)
    σ_hyd =  @. arrh_wet(S_hyd, H_hyd, boltz_k, m.T, m.Ch2o_opx, a, 1f0)

    σ = @. σ_pol +σ_hyd

    return RockphyCond(log10.(σ))
end


mutable struct Dai_Karato2009{F1, F2} <: AbstractCondModel
    T::F1
    Ch2o_opx::F2
end

params_Dai_Karato2009 = (
    A = 10f0 ^2.4f0,
    Aw=10f0 ^ 2.6f0,
    H=147f0,
    Hw=82f0,
    r=62f-2
)

default_params_Dai_Karato2009 =deepcopy(params_Dai_Karato2009)

function forward(m::Dai_Karato2009, p, params = default_params_Zhang2012)
    @unpack A, Aw, H, Hw, r = params
    
    σ_dry =  @. arrh_dry(A, H, gas_R, m.T)
    σ_wet =  @. arrh_wet(Aw, Hw, gas_R_k, m.T, m.Ch2o_opx, 0f0, r)

    σ = @. σ_dry +σ_wet

    return RockphyCond(log10.(σ))
end

mutable struct Yang2011{F1, F2} <: AbstractCondModel
    T::F1
    Ch2o_cpx::F2
end

params_Yang2011 = (
    A = 10f0 ^2.16f0,
    Aw=10f0 ^ 3.56f0,
    H=1.06f0,
    Hw=0.73f0,
    r=113f-2
)

default_params_Yang2011 =deepcopy(params_Yang2011)

function forward(m::Yang2011, p, params = default_params_Yang2011)
    @unpack A, Aw, H, Hw, r = params
    
    σ_dry =  @. arrh_dry(A, H, boltz_k, m.T)
    σ_wet =  @. arrh_wet(Aw, Hw, boltz_k, m.T, m.Ch2o_opx, 0f0, r)

    σ = @. σ_dry +σ_wet

    return RockphyCond(log10.(σ))
end



new_fn(::Nothing) = nothing

for k in subtypes(AbstractCondModel)
    @show k
    s_ = "new_fn(::Type{T}) where {T <: $k} = MT.default_params_$k"
    @show s_
    eval(Meta.parse("new_fn(::Type{T}) where {T <: $k} = MT.default_params_$k"))
end


for k in subtypes(AbstractCondModel)
    @show k
    s_ = "dp_$k = MT.default_params_$k"
    @show s_
    eval(Meta.parse("new_fn(::Type{T}) where {T <: $k} = MT.default_params_$k"))
end

default_params(SEO3)



# =====

abstract type AbstractRockphyModel end
abstract type AbstractRockphyResponse end


# =====


:(mutable struct model1{F}
    T::F
end)


expr1 = :(mutable struct model2{F}
    T::F
end)

expr2 = :(mutable struct model1{F1 <: Union{Distribution, AbstractArray}, F2 <: Union{Distribution, AbstractArray}}
    T::F1
    ch2o::F2
end)


expr1.head
expr1.args
expr1.args[1]
expr1.args[2]
expr1.args[3]

eval(expr1)

model2

expr2
expr2.head
expr2.args
expr2.args[1]
expr2.args[2]

expr2.args[2].head
expr2.args[2].args

expr2.args[2].args[2]
expr2.args[2].args[2].head
expr2.args[2].args[2].args

expr2.args[3]

expr2.args[3].head
expr2.args[3].args

expr2.args[3].args[2]
expr2.args[3].args[2].head
expr2.args[3].args[2].args


# ===

expr_dist = Expr(:curly, :Union, Distribution, AbstractArray)
k_type = subtypes(AbstractCondModel)[1]

fnames = fieldnames(k_type)

expr_type = Expr(:curly, [Symbol("$(k_type)Distribution6"),
    [
        Expr(:<:, Symbol("F$i"), expr_dist) for i in eachindex(fnames)
    ]...
]...)

expr_body = Expr(:block,
    [
        Expr(:(::), Symbol("$(fnames[i])"), Symbol("F$i")) for i in eachindex(fnames)
    ]...
)

expr_3 = Expr(:struct, true, expr_type, expr_body)
eval(expr_3)




docstring = """
$(k_type)Distribution 

Model Distribution for `$k_type`.

## Arguments

Same as `$k_type`
"""

doc_expr = Expr(:macrocall, Symbol("@doc"), nothing, docstring, expr_3)

@eval @doc docstring expr_3

eval(doc_expr)



for k in 8:9
    k_type = subtypes(AbstractCondModel)[k]
    fnames = fieldnames(k_type)

    expr_type = Expr(:curly, [Symbol("$(k_type)Distribution5"),
        [
            Expr(:<:, Symbol("F$i"), expr_dist) for i in eachindex(fnames)
        ]...
    ]...)

    expr_body = Expr(:block,
        [
            Expr(:(::), Symbol("$(fnames[i])"), Symbol("F$i")) for i in eachindex(fnames)
        ]...
    )

    expr_struct2 = Expr(:struct, true, expr_type, expr_body)

    args_block = split(string(Base.Docs.doc(k_type)), "##")[2]

    docstring = """
    $(k_type)Distribution 

    Model Distribution for `$k_type`.

    ## $args_block

    ## Usage
    Refer to the documentation for usage examples.
    """

    # @show quote
    #     @doc $docstring expr_struct2
    # end

    doc_expr = Expr(:macrocall, Symbol("@doc"), nothing, docstring, expr_struct2)

    # Meta.eval(quote
    #         @doc $docstring $expr_struct2
    #     end)

    # Meta.eval(doc_expr)
    # @eval @doc $docstring $expr_struct2
    Meta.eval(Expr(:macrocall, Symbol("@doc"), nothing, docstring, expr_struct))
end



[
        (:("$k")::"F$k") for k in eachindex(fnames)
    ]

for i in eachindex(fnames)
    Expr(:<:, Symbol("F$i"), expr_dist)

end


exp_ = :(Union{Distribution, AbstractArray})

exp_.head
exp_.args


str_ = """
    SEO3Distribution(T)

Model Distribution for `SEO3`[@ref].

## Arguments

  - `T` : Temperature of olivine (in K)

## Usage

Refer to the documentation for usage examples.
"""

split(str_, "##")[2]


expr4 = quote
    sample_type(d::$(Symbol(dstring))) = $mstring
end

eval(expr4)

## ===

macro vectorize_struct(struct_definition)
    # 1. Use MacroTools to easily parse the input expression.
    # This captures the name, parameters, supertype, and fields.
    parts = @capture(struct_definition, struct Name_{Params_} <: SuperType_ Fields__ end)

    # If the capture fails (e.g., wrong input format), do nothing.
    if !parts
        error("Usage: @vectorize_struct struct ... end")
    end

    # 2. Generate the new names for the vectorized version.
    vec_name = Symbol(string(Name), "_vec")
    vec_supertype = Symbol(string(SuperType), "_vec")

    # 3. Create the new parameter definitions.
    # This example makes every parameter a subtype of AbstractArray.
    vec_params = [:($p <: AbstractArray) for p in Params]

    # 4. Generate the docstring for the new struct.
    doc = """
        Auto-generated vectorized version of `$Name`.
        All type parameters have been constrained to `<: AbstractArray`.
        """

    # 5. Build the full expression to be returned.
    # The `quote ... end` block lets us return multiple expressions.
    quote
        # First, make sure the original struct gets defined.
        $(esc(struct_definition))

        # Now, define the new vectorized struct.
        @doc $doc
        struct $(esc(vec_name)){$(esc.(vec_params)...)} <: $(esc(vec_supertype))
            $(esc.(Fields)...)
        end
    end
end

# abstract type AbstractRockphyModelDistribution end
expr = quote
x = subtypes(Real)
for xi in x 
    println(xi)
end
end

expr.head
expr.args

expr.args[4]
expr.args[4].head
expr.args[4].args

macro m1(x)
    # typex = :(typeof($x))
    # Symbol("$(x)Distribution")
    # :(subtypes($x))

    expr = quote
        v = subtypes($x)

        for k in v
            # @show @m2 k

            exprr = :(@get_dist_struct2((Symbol("$k"), fieldnames(k))))
            fnames = fieldnames(k)
            @show exprr
            # Expr(:macrocall, Symbol("@get_dist_struct2"), Symbol("$k"), ("$fnames"))
        end

    end
    esc(expr)
end

macro m5(ktypes, fnames...)
    @show fnames
    # expr = quote
    #     v = fieldnames($x)

    # end
    macro_expr1_ = Expr[]
    fnames2 = esc(:(fieldnames($ktypes)))

    @show fnames2

    for k in eachindex(fnames)
        expr1_ = :(println("this macro => "))

        push!(macro_expr1_, expr1_)
        
    end

    return esc(Expr(:block, macro_expr1_...))

end

# @macroexpand @m1 AbstractCondModel
@m5 UHO2014 :T, :a, :g

@macroexpand @m5 UHO2014 :T :a :g

m = SEO3([100.]);

@m1 m

@macroexpand @m1 m
@macroexpand @m1 SEO3

macro get_dist_struct(k_type, docs1, fnames...)

    expr_dist = Expr(:curly, :Union, Distribution, AbstractArray)
    # expr_dist = :(Union{Distribution, AbstractArray})
    # fnames = :(fieldnames($k_type))

    # @show fnames
    # expr_body = [Expr(:(::), Symbol("$(fnames[i])"), Symbol("F$i")) for i in eachindex(fnames)]
    # @show expr_body

    # struct_name_ = Symbol("$(k_type)Distribution")

    expr_type = Expr(:<:, Expr(:curly, [Symbol("$(k_type)Distribution"),
        [
            Expr(:<:, Symbol("F$i"), expr_dist) for i in eachindex(fnames)
        ]...
    ]...), AbstractRockphyModelDistribution)

    # @show expr_type

    expr_body = Expr(:block,
        [
            Expr(:(::), Symbol("$(fnames[i])"), Symbol("F$i")) for i in eachindex(fnames)
        ]...
    )

    # @show expr_body

    expr_struct = Expr(:struct, true, expr_type, expr_body)

    # @show string(Base.Docs.doc("$k_type"))
    # args_block = split(string(Base.Docs.doc(k_type)), "##")[2]
    docstring = """
    $(k_type)Distribution 

    Model Distribution for `$k_type`.

    ## Arguments

    Same as `$k_type`

    ## Usage
    Refer to the documentation for usage examples.
    """

    doc_expr = Expr(:macrocall, Symbol("@doc"), nothing, docstring, expr_struct)
    # expr_whole = Expr(:macrocall, )
    # expr_whole = :(@doc $doc_expr $expr_struct)

    esc(doc_expr)

end


expr_macro = :(@get_dist_struct(SEO3, :T))


@macroexpand @get_dist_struct SEO3 :T
@get_dist_struct SEO3 :T
# @macroexpand @get_dist_struct SEO3

k_type1 = [Sifre2014, Gaillard2008]
fnames1 = fieldnames.(k_type1)

@get_dist_struct2(k_type1, fnames1...)
@macroexpand @get_dist_struct2(k_type1, fnames1...)
@macroexpand @get_dist_struct2 k_type1 fnames1...

for i in 1:2
    k_ = k_type1[i]
    f_ = fnames1[i]
    expr_m =quote
        @get_dist_struct2($k_, $(f_...))
    end
    @show expr_m
end


exp_1 = Expr(:block, expr_block[1])

ee_ = @eval exp_1

ex_1 = expr_block[1]

@eval ex_1

# ====



expr = :(struct Foo
    a::Int = 1         # specified default
    b::String          # required keyword
end)

@macroexpand @kwdef struct Foo
    a::Int = 1         # specified default
    b::String          # required keyword
end

# expr = quote
    #     v = fieldnames($x)

# end

macro distdef(ex)
    expr_ = quote
        # ex

        struct_name_ = "$ex.args[2]"*"Distribution"
        expr_type = Expr(:<:, Expr(:curly, [Symbol(struct_name_),
      [
          Expr(:<:, Symbol("F$i"), expr_dist) for i in eachindex(fnames)
      ]...
  ]...), AbstractRockphyModelDistribution)
    end
end

expr_dist = Expr(:curly, :Union, Distribution, AbstractArray)

function make_struct_expr(k_type, fnames)
    struct_name_ = Symbol("$(k_type)Distribution")

  expr_type = Expr(:<:, Expr(:curly, [Symbol("$(k_type)Distribution"),
      [
          Expr(:<:, Symbol("F$i"), expr_dist) for i in eachindex(fnames)
      ]...
  ]...), AbstractRockphyModelDistribution)

  expr_body = Expr(:block,
      [
          Expr(:(::), Symbol("$(fnames[i])"), Symbol("F$i")) for i in eachindex(fnames)
      ]...
  )

  expr_struct_ = Expr(:struct, true, expr_type, expr_body)

  args_block = split(string(Base.Docs.doc(k_type)), "##")[2]

    docstring = """
    $(k_type)Distribution 

    Model Distribution for `$k_type`.

    ## $args_block

    ## Usage
    Refer to the documentation for usage examples.
    """

    doc_expr = Expr(:macrocall, Symbol("@doc"), nothing, docstring, expr_struct_)

end

function make_docstring(k_type)
    
    args_block = split(string(Base.Docs.doc(k_type)), "##")[2]

    docstring = """
    $(k_type)Distribution 

    Model Distribution for `$k_type`.

    ## $args_block

    ## Usage
    Refer to the documentation for usage examples.
    """

    doc_expr = Expr(:macrocall, Symbol("@doc"), nothing, docstring, expr_struct_)

end

make_struct_expr(SEO3, fieldnames(SEO3))

make_docstring(SEO3)

macro make_dist(k_type, fnames)
    @show k_type
    @show fnames
    expr = :(make_struct_expr($k_type, $fnames))
    esc(@eval $expr)
end

k_type_ = SEO3
fnames_ = fieldnames(SEO3)

@make_dist k_type_ fnames_

@macroexpand @make_dist k_type_ fnames_

exp1 = :(
    global mutable struct c1
        T
    end
)

exp2_ = :(
    mutable struct c1
        T
    end
)


exp1.head
exp1.args

Expr(:global, exp2_)
