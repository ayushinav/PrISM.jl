import SubsurfaceCore: plot_model!
function SubsurfaceCore.plot_model!(
        ax, model::AbstractGeophyModel{T1, T2}; half_space_thickness=nothing,
        kwargs...) where {T1 <: AbstractVector, T2 <: GridapType}
    @show "HERE"

    zgrid_ = GridapGmsh.get_node_coordinates(Triangulation(model.h))
    zgrid = [z.data[1] for z in zgrid_]

    stairs!(ax, [mgrid..., mgrid[end]], zgrid; kwargs...)
    ax.yreversed = true
end
