using DrWatson
@quickactivate
include(srcdir("Plot_utilities.jl"))

using Bonito, Observables
using TidierDB, WGLMakie, JLD2, DataFrames, StatsBase, Statistics, CSV
using Memoize
import WGLMakie as W

# ==========================================
# 1. Data and Helper Definitions
# ==========================================

@load "data/dipterans/D_melanogaster_GENE_DB.jld2" GENE_DB gene_table
con = DBInterface.connect(DuckDB.DB("data/ERC.duckdb"))
dros_table = db_table(con, "D_melanogaster")

ERC_total = @chain dros_table begin
    @select(fERC)
    @collect
end
ERC_total = ERC_total.fERC

const JLD2_CACHE_DIR = joinpath(@__DIR__, "data", "jld2_cache")
mkpath(JLD2_CACHE_DIR)

function disk_cache(f::Function, category::String, key::String)
    cat_dir = joinpath(JLD2_CACHE_DIR, category)
    mkpath(cat_dir)
    file_path = joinpath(cat_dir, "$(key).jld2")
    if isfile(file_path)
        return JLD2.load_object(file_path)
    else
        result = f()
        JLD2.save_object(file_path, result)
        return result
    end
end

function make_sortable_th(label::String; style::String="")
    base_style = "padding: 8px; background-color: #f2f2f2; border: 1px solid #ddd; cursor: pointer; user-select: none; " * style
    return DOM.th(
        label,
        style = base_style,
        onclick = js"""function(event) {
            let th = event.currentTarget;
            let table = th.closest('table');
            let tbody = table.querySelector('tbody');
            if (!tbody) return;

            let rows = Array.from(tbody.querySelectorAll('tr'));
            let index = Array.from(th.parentNode.children).indexOf(th);
            let isAsc = th.dataset.sortDir !== 'asc';

            // Clear sort state and arrow indicators on sibling headers
            th.parentNode.querySelectorAll('th').forEach(h => {
                h.dataset.sortDir = '';
                let arrow = h.querySelector('.sort-arrow');
                if (arrow) arrow.remove();
            });

            th.dataset.sortDir = isAsc ? 'asc' : 'desc';

            // Append directional arrow indicator
            let arrow = document.createElement('span');
            arrow.className = 'sort-arrow';
            arrow.textContent = isAsc ? ' ▲' : ' ▼';
            arrow.style.fontSize = '11px';
            arrow.style.marginLeft = '4px';
            arrow.style.color = '#333';
            th.appendChild(arrow);

            // Sort rows based on data type (numeric vs string)
            rows.sort((a, b) => {
                let valA = a.children[index] ? a.children[index].innerText.trim() : '';
                let valB = b.children[index] ? b.children[index].innerText.trim() : '';

                let numA = parseFloat(valA);
                let numB = parseFloat(valB);

                let cmp;
                if (!isNaN(numA) && !isNaN(numB) && valA !== '' && valB !== '') {
                    cmp = numA - numB;
                } else {
                    cmp = valA.localeCompare(valB, undefined, { numeric: true, sensitivity: 'base' });
                }
                return isAsc ? cmp : -cmp;
            });

            // Re-append sorted rows to the table body
            rows.forEach(r => tbody.appendChild(r));
        }"""
    )
end

function get_precomputed_stats(gene)
    idx = findfirst(gene_table.gene .== gene)
    return Dict(
        "Name/Symbol" => gene_table.name[idx],
        "FlybaseID" => gene_table.stable_id[idx],
        "Position" => gene_table.chr[idx] * ":" * string(gene_table.pos_1[idx]) * "-" * string(gene_table.pos_2[idx]), 
        "Description" => gene_table.description[idx],
        "mean fERC" => gene_table.mean_fERC[idx],
        "RBH in data" => gene_table.num_rbh[idx],
        "Cluster in network" => gene_table.community[idx]
    )
end

# Memoized: Pure DataFrame computation 
@memoize function get_interactions(gene)
    return disk_cache("interactions", gene) do
        ERC_table = @eval @chain dros_table begin
            @filter(i == $gene || j == $gene)
            @select(i, j, r, fERC, pval,rfERC,rpval)
            @collect
        end
        swaps = findall(ERC_table.j .== gene)
        ERC_table.j[swaps] = ERC_table.i[swaps]
        ERC_table.i .= gene
        ERC_table.PartnerName = [coalesce(gene_table.name[findfirst(gene_table.gene .== x)], "") for x in ERC_table.j]
        return sort(ERC_table,:fERC,rev=true)
    end
end

# Memoized: Heavy DuckDB queries and vector transformations ONLY, cached for later searches
@memoize function get_plot_data(gene)
    return disk_cache("plot_data", gene) do
        chr_order = ["2L", "2R", "3L", "3R", "4", "X", "Y", "mitochondrion_genome"]
        chr_lengths = [23312347, 25254477, 28012551, 32067996, 1271129, 23107330, 3666928]

        ERC_data = @eval @chain dros_table begin
            @filter(i == $gene || j == $gene)
            @select(i, j, fERC, same_chr, dist, chr_i, chr_j, pos_i_1, pos_i_2, pos_j_1, pos_j_2, pval, null)
            @collect
        end

        idx = findall(ERC_data.i .== gene)
        partners = copy(ERC_data.i)
        partners[idx] = ERC_data.j[idx]
        
        partner_names = [coalesce(gene_table.name[findfirst(gene_table.gene .== x)], " ") for x in partners]
        chrs = copy(ERC_data.chr_i)
        chrs[idx] = ERC_data.chr_j[idx]

        pos1 = copy(ERC_data.pos_i_1)
        pos2 = copy(ERC_data.pos_i_2)
        pos1[idx] = ERC_data.pos_j_1[idx]
        pos2[idx] = ERC_data.pos_j_2[idx]

        pos = (pos1 .+ pos2) ./ 2
        pos_converted = chr_order_pos(pos, chrs, chr_order, chr_lengths)
        mean_fERC = [gene_table.mean_fERC[findfirst(gene_table.gene .== x)] for x in partners]
        ratio = ERC_data.fERC .- mean_fERC

        return (
            ERC_data = ERC_data,
            partners = partners,
            partner_names = partner_names,
            pos_converted = pos_converted,
            mean_fERC = mean_fERC,
            ratio = ratio
        )
    end
end

#Can't memoize since non-static objects
function makie_plot(gene)
    data = get_plot_data(gene)

    chr_order = ["2L", "2R", "3L", "3R", "4", "X", "Y", "mitochondrion_genome"]
    chr_lengths = [23312347, 25254477, 28012551, 32067996, 1271129, 23107330, 3666928]
    chr_lines = [sum(chr_lengths[1:x]) for x in 1:7]
    chr_label_pos = vcat(chr_lengths[1]/2, [chr_lengths[i]/2 + chr_lines[i-1] for i in 2:7])
    chr_ticks = (chr_label_pos, chr_order[1:7])

    fig = W.Figure()
    
    ax1 = Axis(fig[1, 1], xlabel="fERC", ylabel="Normalized Frequency", xgridvisible=false, ygridvisible=false)
    density!(ax1, data.ERC_data.fERC, strokecolor=:red, color=(:red, 0.5), label="Gene")
    density!(ax1, ERC_total, strokecolor=:blue, color=(:blue, 0.5), label="Overall")
    xlims!(ax1, (-5, 15))
    vlines!(ax1, quantile(data.ERC_data.fERC, [0.5]), color=:red)
    vlines!(ax1, quantile(ERC_total, [0.5]), color=:blue)

    ax2 = Axis(fig[1, 2], xlabel=L"\overline{fERC}", ylabel=L"fERC_{gene}", xgridvisible=false, ygridvisible=false)
    sc2 = scatter!(ax2, data.mean_fERC, data.ERC_data.fERC, color=data.ratio, markersize=3, colormap=:vik,
                   inspector_label=(self, index, position) -> "Gene_ID: " * data.partners[index] * "\n Name: " * data.partner_names[index])
    lines!(ax2, range(-2, 6), range(-2, 6), linestyle=:dash, color=:black)
    xlims!(ax2, [-2, 5])
    Colorbar(fig[1, 2][1, 2], sc2, label=L"fERC_{gene}-\overline{fERC}")

    ax3 = Axis(fig[2, :], xlabel="Position", ylabel=L"fERC_{Sxl}", yminorgridvisible=false, xticks=chr_ticks,
               xticksize=0, xminorticks=chr_lines, xminorticksvisible=true, xminorgridvisible=true,
               xminorgridcolor=RGBAf(0, 0, 0, 0.25), xgridcolor=RGBAf(0, 0, 0, 0),
               limits=(0, maximum(chr_lines), nothing, nothing))

    sc3 = scatter!(ax3, data.pos_converted, data.ERC_data.fERC, color=data.ratio, markersize=3, colormap=:vik,
                   inspector_label=(self, index, position) -> "Gene_ID: " * data.partners[index] * "\n Name: " * data.partner_names[index])
    
    DataInspector(fig)
    return fig
end

# Memoized: GO CSV lookup. Cost - can't sort columns! 
@memoize function go_view(gene)
    file_path = "data/dipterans/D_melanogaster/gene_enrichment/" * gene * ".csv"
    if isfile(file_path)
        go_res = CSV.read(file_path, DataFrame)
        rows = [DOM.tr(
            DOM.td(x.GO_ID, style="padding: 8px; border: 1px solid #ddd;"),
            DOM.td(string(x.pval_BH), style="padding: 8px; border: 1px solid #ddd;"),
            DOM.td(string(x.exp_fERC), style="padding: 8px; border: 1px solid #ddd;"),
            DOM.td(string(x.obs_fERC), style="padding: 8px; border: 1px solid #ddd;"),
            DOM.td(x.description, style="padding: 8px; border: 1px solid #ddd;")
        ) for x in eachrow(go_res)]
        
        header = DOM.tr(
            make_sortable_th("GO_ID"), 
            make_sortable_th("adjusted pval"),
            make_sortable_th("expected fERC"),
            make_sortable_th("observed fERC"),
            make_sortable_th("description")
        )
        return DOM.div(
            DOM.h3("GO Depletion/Enrichment"), 
            DOM.table(
                DOM.thead(header), 
                DOM.tbody(rows...), 
                style="width: 100%; border-collapse: collapse; text-align: left;"
            ),
            style="overflow-y:auto; height:600px;"
        )
    else
        return DOM.div(DOM.h3("No GO Enrichment/Depletion File"))
    end
end

function stats_view(gene)
    stats = get_precomputed_stats(gene)
    items = [DOM.li(DOM.strong("$k: "), "$v") for (k, v) in stats]
    return DOM.div(DOM.h3("Summary Statistics"), DOM.ul(items...))
end

function interactions_view(gene) 
    interactions = get_interactions(gene)
    if size(interactions)[1] > 0
        rows = [
            DOM.tr(
                DOM.td(x.j, style="padding: 8px; border: 1px solid #ddd;"),
                DOM.td(x.PartnerName, style="padding: 8px; border: 1px solid #ddd;"),
                DOM.td(string(x.r), style="padding: 8px; border: 1px solid #ddd;"),
                DOM.td(string(x.fERC), style="padding: 8px; border: 1px solid #ddd;"),
                DOM.td(string(x.pval), style="padding: 8px; border: 1px solid #ddd;"),
                DOM.td(string(x.rfERC), style="padding: 8px; border: 1px solid #ddd;"),
                DOM.td(string(x.rpval), style="padding: 8px; border: 1px solid #ddd;")
            ) for x in eachrow(interactions)
        ]

        header = DOM.tr(
            make_sortable_th("Partner ID", style="position: sticky; top: 0; z-index: 1;"),
            make_sortable_th("Partner Name", style="position: sticky; top: 0; z-index: 1;"),
            make_sortable_th("Pearson Correlation", style="position: sticky; top: 0; z-index: 1;"),
            make_sortable_th("fERC", style="position: sticky; top: 0; z-index: 1;"),
            make_sortable_th("p-val (correlation)", style="position: sticky; top: 0; z-index: 1;"),
            make_sortable_th("rfERC", style="position: sticky; top: 0; z-index: 1;"),
            make_sortable_th("p-val (residual)", style="position: sticky; top: 0; z-index: 1;")
        )

        html_table = DOM.table(
            DOM.thead(header), 
            DOM.tbody(rows...), 
            style="width: 100%; border-collapse: collapse; text-align: left;"
        )

        # Checkbox toggle control
        filter_toggle = DOM.div(
            DOM.label(
                DOM.input(
                    type = "checkbox",
                    onchange = js"""function(event) {
                        let isChecked = event.target.checked;
                        let container = event.target.closest('.interactions-container');
                        let rows = container.querySelectorAll('tbody tr');
                        
                        rows.forEach(r => {
                            // Column index 4 corresponds to p-val, adjusted threshold is 1e-7
                            let pvalCell = r.children[4];
                            if (pvalCell) {
                                let pval = parseFloat(pvalCell.innerText.trim());
                                if (isChecked && !isNaN(pval) && pval >= 1e-7) {
                                    r.style.display = 'none';
                                } else {
                                    r.style.display = '';
                                }
                            }
                        });
                    }""",
                    style = "margin-right: 8px; cursor: pointer;"
                ),
                "Show only significant interactions (p_adjusted < 0.05)",
                style = "cursor: pointer; user-select: none; font-size: 14px; font-weight: 500;"
            ),
            style = "margin-bottom: 12px; display: flex; align-items: center;"
        )

        return DOM.div(
            DOM.h3("Interactions"),
            filter_toggle,
            DOM.div(
                html_table,
                style="overflow-y:auto; height:300px;"
            ),
            class = "interactions-container"
        )
    else 
        return DOM.div(DOM.h3("No interaction data"))
    end
end

function figure_view(gene)
    # 1. Create the Observable holding the loading UI
    plot_container = Observable{Any}(
        DOM.div(
            DOM.style("""
                @keyframes spin {
                    0% { transform: rotate(0deg); }
                    100% { transform: rotate(360deg); }
                }
                .spinner {
                    border: 4px solid #f3f3f3;
                    border-top: 4px solid #3498db;
                    border-radius: 50%;
                    width: 36px;
                    height: 36px;
                    animation: spin 1s linear infinite;
                    margin: 0 auto 12px auto;
                }
            """),
            DOM.div(
                DOM.div(class="spinner"),
                DOM.p("Generating plot for ", DOM.strong(gene), "... Please wait.", 
                      style="color: #555; font-size: 15px; margin: 0;"),
                style="text-align: center; padding: 40px; background: #fafafa; border: 1px dashed #ccc; border-radius: 8px; margin: 20px 0;"
            )
        )
    )

    # 2. Asynchronously compute plot in background
    @async begin
        try
            fig = makie_plot(gene)
            plot_container[] = DOM.div(fig)
        catch e
            @error "Failed to generate plot for $gene" exception=(e, catch_backtrace())
            plot_container[] = DOM.div(
                "An error occurred while generating the plot: $e",
                style="color: #d9534f; padding: 20px; background: #fdf2f2; border-radius: 6px;"
            )
        end
    end

    # 3. Return a DOM Node wrapping the Observable
    return DOM.div(plot_container)
end

# ==========================================
# 2. App Definition
# ==========================================

dashboard = App() do session
    search_query = Observable("")
    selected_gene = Observable("NP_001259303.1")
    active_tab = Observable("Overview") 
    
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
            contains(lowercase(coalesce(g.symbol, "")), q) || 
            contains(lowercase(coalesce(g.id, "")), q) ||
            contains(lowercase(coalesce(g.alias, "")), q)
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
                    $(search_query).notify(''); 
                    $(active_tab).notify('Overview');
                }"""
            )
        end
        
        return DOM.div(
            items..., 
            style="border: 1px solid #ccc; border-top: none; border-radius: 0 0 4px 4px; max-height: 400px; overflow-y: auto; box-shadow: 0 4px 6px rgba(0,0,0,0.1);"
        )
    end
    
    header_view = map(selected_gene) do gene
        DOM.h2("Currently Viewing: $gene", style="color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px;")
    end

    dynamic_view = map(active_tab, selected_gene) do tab, gene
        if tab == "Overview"
            return stats_view(gene)
        elseif tab == "Plot"
            return DOM.div(DOM.h3("Plot View"),DOM.tbody("This is an interactive figure. Hover over points to identify genes.") ,figure_view(gene))
        elseif tab == "Interactions"
            return interactions_view(gene)
        elseif tab == "GO"
            return go_view(gene)
        else 
            return stats_view(gene)
        end
    end

    tab_bar = DOM.div(
        DOM.button("Overview", onclick=js"function(event) { $(active_tab).notify('Overview'); }"),
        DOM.button("Plot", onclick=js"function(event) { $(active_tab).notify('Plot'); }"),
        DOM.button("Interactions", onclick=js"function(event) { $(active_tab).notify('Interactions'); }"),
        DOM.button("GO", onclick=js"function(event) { $(active_tab).notify('GO'); }"),
        style="display: flex; gap: 10px; margin-bottom: 20px;"
    )

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

server = Bonito.Server(dashboard, "0.0.0.0", 8080)