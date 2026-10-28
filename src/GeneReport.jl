using DrWatson
@quickactivate
include(srcdir("Plot_utilities.jl"))
using Bonito, Observables

using TidierDB, WGLMakie, JLD2, DataFrames, StatsBase, Statistics, CSV
import WGLMakie as W
# ==========================================
# 1. Reading in data and function definitions.
# ==========================================


@load "data/dipterans/D_melanogaster_GENE_DB.jld2" GENE_DB gene_table
con = DBInterface.connect(DuckDB.DB("data/ERC.duckdb"))
dros_table = db_table(con,"D_melanogaster")
ERC_total = @chain dros_table begin
    @select(fERC)
    @collect
end
ERC_total = ERC_total.fERC

function get_precomputed_stats(gene)
    idx = findfirst(gene_table.gene .== gene)
    return(Dict("Name/Symbol"=>gene_table.name[idx],
                "FlybaseID"=>gene_table.stable_id[idx],
                "Position" => gene_table.chr[idx]*":"*string(gene_table.pos_1[idx])*"-"*string(gene_table.pos_2[idx]), 
                "Description" =>gene_table.description[idx],
                "mean fERC"=>gene_table.mean_fERC[idx],
                "RBH in data"=>gene_table.num_rbh[idx],
                "Cluster in network"=>gene_table.community[idx]))
end
function get_significant_interactions(gene)
    idx = findfirst(gene_table.gene .==gene)
    partners = gene_table.partners[idx]
    partner_names = gene_table.name[[findfirst(gene_table.gene .== x) for x in partners]]
    return(DataFrame(PartnerID=partners,PartnerName=partner_names,fERC=gene_table.partner_fERC[idx]))
end

function makie_plot(gene)
    chr_order = ["2L","2R","3L","3R","4","X","Y","mitochondrion_genome"]
    chr_lengths = [23312347,25254477,28012551,32067996,1271129,23107330,3666928]
    chr_lines = [sum(chr_lengths[1:x]) for x in 1:7]
    chr_label_pos = vcat(chr_lengths[1]/2,[chr_lengths[i]/2+chr_lines[i-1] for i in 2:7])
    chr_ticks = (chr_label_pos,chr_order[1:7])

    ERC_data = @eval @chain dros_table begin
               @filter(i == $gene || j == $gene)
               @select(i,j,fERC,same_chr,dist,chr_i,chr_j,pos_i_1,pos_i_2,pos_j_1,pos_j_2,pval,null)
               @collect
    end
    fig = W.Figure()
    ax1 = Axis(fig[1,1],
    xlabel="fERC",
    ylabel="Normalized Frequency",
    xgridvisible=false,
    ygridvisible=false)
    density!(ax1,ERC_data.fERC,strokecolor=:red,color=(:red,0.5),label="Gene")
    density!(ax1,ERC_total,strokecolor=:blue,color=(:blue,0.5),label="Overall")
    xlims!(ax1,(-5,15))
    vlines!(ax1,quantile(ERC_data.fERC,[0.5]),color=:red)
    vlines!(ax1,quantile(ERC_total,[0.5]),color=:blue)

    
    idx = findall(ERC_data.i .== gene)
    partners = ERC_data.i
    partners[idx] = ERC_data.j[idx]
    partner_names = [coalesce(gene_table.name[findfirst(gene_table.gene .== x)]," ") for x in partners]
    chrs = ERC_data.chr_i
    chrs[idx] = ERC_data.chr_j[idx]

    pos1 = ERC_data.pos_i_1
    pos2 = ERC_data.pos_i_2
    pos1[idx] = ERC_data.pos_j_1[idx]
    pos2[idx] = ERC_data.pos_j_2[idx]

    pos = (pos1+pos2)/2

    pos_converted = chr_order_pos(pos,chrs,chr_order,chr_lengths)
    mean_fERC = [gene_table.mean_fERC[findfirst(gene_table.gene .== x)] for x in partners]
    ERC_data.ratio = ERC_data.fERC .- mean_fERC

    
    ax2 = Axis(fig[1,2],
    xlabel=L"\overline{fERC}",
    ylabel=L"fERC_{gene}",
    xgridvisible=false,
    ygridvisible=false)
    sc2= scatter!(ax2,mean_fERC,ERC_data.fERC,color=ERC_data.ratio,markersize=3,colormap=:vik,inspector_label=(self,index,position) -> "Gene_ID: "*partners[index]*"\n Name: "*partner_names[index])
    lines!(ax2,range(-2,6),range(-2,6),linestyle=:dash,color=:black)
    xlims!(ax2,[-2,5])
    Colorbar(fig[1,2][1,2],sc2,label=L"fERC_{gene}-\overline{fERC}")


    ax3 = Axis(fig[2,:],
    xlabel="Position",
    ylabel=L"fERC_{Sxl}",
    yminorgridvisible=false,
    xticks=chr_ticks,
    xticksize=0,
    xminorticks=chr_lines,
    xminorticksvisible=true,
    xminorgridvisible=true,
    xminorgridcolor=RGBAf(0,0,0,0.25),
    xgridcolor=RGBAf(0,0,0,0),
    limits =(0,maximum(chr_lines),nothing,nothing))

    sc3 = scatter!(ax3,pos_converted,ERC_data.fERC,color=ERC_data.ratio,markersize=3,colormap=:vik,inspector_label=(self,index,position) -> "Gene_ID: "*partners[index]*"\n Name: "*partner_names[index])
    DataInspector(fig)
    return(fig)
end

# ==========================================
# 2. App Definition
# ==========================================

dashboard = App() do session
    search_query = Observable("")
    selected_gene = Observable("NP_001259303.1")
    active_tab = Observable("Overview") 
    
    # ------------------------------------------
    # Search Box & Autocomplete UI
    # ------------------------------------------
    
    search_box = DOM.input(
        id = "gene-search",
        value = search_query,
        type = "text",
        placeholder = "Search by gene symbol, ID, or alias (e.g., p53)...",
        oninput = js"function(event) { $(search_query).notify(event.target.value); }",
        style = "width: 100%; padding: 12px; font-size: 16px; border: 1px solid #ccc; border-radius: 4px; box-sizing: border-box;"
    )

    search_results = map(search_query) do query
        if length(query) < 2
            return DOM.div("Type at least 2 characters to search...", style="color: #888; padding: 10px;")
        end
        
        q = lowercase(query)
        
        matches = filter(GENE_DB) do g
            contains(lowercase(coalesce(g.symbol,"")), q) || 
            contains(lowercase(coalesce(g.id,"")), q) ||
            contains(lowercase(coalesce(g.alias,"")),q)
            #The below is once I annotate more aliases per gene, for now just fbid
            #any(a -> contains(lowercase(a), q), g.alias)
        end
        
        if isempty(matches)
            return DOM.div("No matching genes found.", style="color: #d9534f; padding: 10px;")
        end
        
        top_matches = first(matches, 10)
        
        items = map(top_matches) do g
            DOM.div(
                DOM.div(DOM.strong(g.symbol), " (", g.id, ")"),
                DOM.div("Aliases: ", join(g.alias, ", "), style="font-size: 12px; color: #666;"),
                style = "padding: 10px; border-bottom: 1px solid #eee; cursor: pointer; background: #fff; transition: background 0.2s;",
                
                onclick = js"""function(event) {
                    $(selected_gene).notify($(g.id));
                    $(search_query).notify(''); // <-- This will now reactively clear the input box via Bonito
                    $(active_tab).notify('Overview')
                }"""
            )
        end
        
        return DOM.div(
            items..., 
            style="border: 1px solid #ccc; border-top: none; border-radius: 0 0 4px 4px; max-height: 400px; overflow-y: auto; box-shadow: 0 4px 6px rgba(0,0,0,0.1);"
        )
    end
    
    # ------------------------------------------
    # Downstream Reactive Components
    # ------------------------------------------

    header_view = map(selected_gene) do gene
        DOM.h2("Currently Viewing: $gene", style="color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px;")
    end

    stats_view = map(selected_gene) do gene
        stats = get_precomputed_stats(gene)
        items = [DOM.li(DOM.strong("$k: "), "$v") for (k, v) in stats]
        return DOM.div(DOM.h3("Summary Statistics"), DOM.ul(items...))
    end
    
    interactions_view = map(selected_gene) do gene
        interactions = get_significant_interactions(gene)
        if size(interactions)[2] > 0
            rows = [
            DOM.tr(
                DOM.td(x.PartnerID, style="padding: 8px; border: 1px solid #ddd;"),
                DOM.td(x.PartnerName, style="padding: 8px; border: 1px solid #ddd;"),
                DOM.td(string(x.fERC), style="padding: 8px; border: 1px solid #ddd;")
            ) for x in eachrow(interactions)
            ]
            header = DOM.tr(
                DOM.th("Partner ID", style="padding: 8px; background-color: #f2f2f2; border: 1px solid #ddd;position: sticky; top: 0; z-index: 1;"),
                DOM.th("Partner Name", style="padding: 8px; background-color: #f2f2f2; border: 1px solid #ddd;position: sticky; top: 0; z-index: 1;"),
                DOM.th("fERC", style="padding: 8px; background-color: #f2f2f2; border: 1px solid #ddd;position: sticky; top: 0; z-index: 1;")
            )
    
            html_table = DOM.table(
                DOM.thead(header), 
                DOM.tbody(rows...), 
                style="width: 100%; border-collapse: collapse; text-align: left;"
            )
            return DOM.div(
            DOM.h3("Significant Interactions",style="position: sticky; top: 0; z-index: 1;"), 
            html_table,
            style="overflow-y:auto; height:300px;"
            )
        else 
            return DOM.div(DOM.h3("No significant interactions"))
        end
    end

    
    figure_view = map(selected_gene) do gene
        return DOM.div(makie_plot(gene))   
    end
    
    go_plot_view = map(selected_gene) do gene
        if isfile("data/dipterans/D_melanogaster/gene_enrichment/"*gene*".csv")
            go_res = CSV.read("data/dipterans/D_melanogaster/gene_enrichment/"*gene*".csv",DataFrame)
            rows = [DOM.tr(DOM.td(x.GO_ID, style="padding: 8px; border: 1px solid #ddd;"),
                       DOM.td(string(x.pval_BH), style="padding: 8px; border: 1px solid #ddd;"),
                       DOM.td(string(x.exp_fERC), style="padding: 8px; border: 1px solid #ddd;"),
                       DOM.td(string(x.obs_fERC), style="padding: 8px; border: 1px solid #ddd;"),
                       DOM.td(x.description, style="padding: 8px; border: 1px solid #ddd;")) for x in eachrow(go_res)]
            header = DOM.tr(DOM.th("GO_ID", style="padding: 8px; background-color: #f2f2f2; border: 1px solid #ddd;"), 
                        DOM.th("adjusted pval", style="padding: 8px; background-color: #f2f2f2; border: 1px solid #ddd;"),
                        DOM.th("expected fERC", style="padding: 8px; background-color: #f2f2f2; border: 1px solid #ddd;"),
                        DOM.th("observed fERC", style="padding: 8px; background-color: #f2f2f2; border: 1px solid #ddd;"),
                        DOM.th("description", style="padding: 8px; background-color: #f2f2f2; border: 1px solid #ddd;")
                        )
            return DOM.div(DOM.h3("GO Depletion/Enrichment"), DOM.table(DOM.thead(header), DOM.tbody(rows...), style="width: 100%; border-collapse: collapse; text-align: left;"),style="overflow-y:auto;height:600px;")
        else
            return DOM.div(DOM.h3("No GO Enrichment/Depletion"))
        end
    end
    
    dynamic_view = map(active_tab) do tab
        if tab == "Overview"
            return stats_view
        elseif tab == "Plot"
            return figure_view
        elseif tab == "Interactions"
            return interactions_view
        elseif tab == "GO"
            return go_plot_view
        else 
            return stats_view
        end
    end

    tab_bar = DOM.div(
        DOM.button("Overview", onclick=onclick = js"""function(event) {
                    $(active_tab).notify("Overview");
                }"""),
        DOM.button("Plot", onclick=js"""function(event) {
                    $(active_tab).notify("Plot");
                }"""),
        DOM.button("Interactions",onclick=js"""function(event) {
                    $(active_tab).notify("Interactions");
                }"""),
        DOM.button("GO", onclick=js"""function(event) {
                    $(active_tab).notify("GO");
                }"""),
        style="display: flex; gap: 10px; margin-bottom: 20px;"
    )

    # ------------------------------------------
    # Layout Assembly
    # ------------------------------------------
    return DOM.div(
        style = "font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; max-width: 1200px; margin: 0 auto; padding: 20px;",
        
        DOM.h1("Co-evolutionary Analytics Dashboard"),
        DOM.p("Search across 12,000+ genes by symbol or Ensembl ID."),
        
        DOM.div(
            search_box,
            search_results,
            style="position: relative; margin-bottom: 40px; max-width: 600px; z-index: 100;"
        ),
        
        header_view,

        tab_bar,

        dynamic_view
    )
end

# ==========================================
# 3. Serve the App
# ==========================================
#server = Bonito.Server(dashboard, "127.0.0.1", 8080)
#println("Dashboard is running! Open http://127.0.0.1:8080 in your browser.")