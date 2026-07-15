using Bonito, Observables
using WGLMakie

using DataFrames, JLD2, CSV
# ==========================================
# 1. Reading in data and making fetch functions.
# ==========================================


@load "data/dipterans/D_melanogaster_GENE_DB.jld2" GENE_DB gene_table

function get_precomputed_stats(gene)
    idx = findfirst(gene_table.gene .== gene)
    return(Dict("Name/Symbol"=>gene_table.name[idx],
                "FlybaseID"=>gene_table.stable_id[idx],
                "Chromosome" => gene_table.chr[idx], 
                "Description" =>gene_table.description[idx],
                "mean fERC"=>gene_table.mean_fERC[idx],
                "RBH in data"=>gene_table.num_rbh[idx],
                "Cluster in network"=>gene_table.community[idx]))
end
function get_significant_interactions(gene)
    idx = findfirst(gene_table.gene .==gene)
    return(DataFrame(partner=gene_table.partners[idx],ferc=gene_table.partner_fERC[idx]))
end

# ==========================================
# 2. App Definition
# ==========================================

dashboard = App() do session
    search_query = Observable("")
    selected_gene = Observable("NP_001259303.1") 
    
    # ------------------------------------------
    # Search Box & Autocomplete UI
    # ------------------------------------------
    
    # FIX 1: Use `event.target.value` instead of `this.value`
    # We also assign an explicit ID so we can clear it easily later.
    search_box = DOM.input(
        id = "gene-search",
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
                
                # FIX 2: Explicitly pass function(event), target the specific input ID, and notify
                onclick = js"""function(event) {
                    $(selected_gene).notify($(g.id));
                    document.getElementById('gene-search').value = '';
                    $(search_query).notify('');
                }""",
                
                onmouseover = js"this.style.background='#f0f8ff'",
                onmouseout = js"this.style.background='#fff'"
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
        rows = [DOM.tr(DOM.td(g.partner, style="padding: 8px; border: 1px solid #ddd;"), DOM.td(string(g.ferc), style="padding: 8px; border: 1px solid #ddd;")) for g in eachrow(interactions)]
        header = DOM.tr(DOM.th("Interacting Gene", style="padding: 8px; background-color: #f2f2f2; border: 1px solid #ddd;"), DOM.th("fERC", style="padding: 8px; background-color: #f2f2f2; border: 1px solid #ddd;"))
        return DOM.div(DOM.h3("Significant Interactions"), DOM.table(DOM.thead(header), DOM.tbody(rows...), style="width: 100%; border-collapse: collapse; text-align: left;"))
    end
    
    figure_view = map(selected_gene) do gene
        img_url = "https://via.placeholder.com/500x350.png?text=$gene+Pre-generated+Figure"
        return DOM.div(DOM.h3("Structural Figure"), DOM.img(src=img_url, style="max-width: 100%; border-radius: 8px;"))
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
            return DOM.div(DOM.h3("GO Depletion/Enrichment"), DOM.table(DOM.thead(header), DOM.tbody(rows...), style="width: 100%; border-collapse: collapse; text-align: left;"))
        else
            return DOM.div(DOM.h3("No GO Enrichment/Depletion"))
        end
    end
    
    # ------------------------------------------
    # Layout Assembly
    # ------------------------------------------
    return DOM.div(
        style = "font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; max-width: 1200px; margin: 0 auto; padding: 20px;",
        
        DOM.h1("Gene Analytics Dashboard"),
        DOM.p("Search across 12,000+ genes by symbol, Ensembl ID, or known aliases."),
        
        DOM.div(
            search_box,
            search_results,
            style="position: relative; margin-bottom: 40px; max-width: 600px; z-index: 100;"
        ),
        
        header_view,
        
        DOM.div(
            DOM.div(stats_view, DOM.br(), interactions_view, style="flex: 1; padding-right: 20px;"),
            DOM.div(figure_view, DOM.br(), go_plot_view, style="flex: 1; display: flex; flex-direction: column; gap: 20px;"),
            style="display: flex; flex-direction: row; flex-wrap: wrap;"
        )
    )
end

# ==========================================
# 3. Serve the App
# ==========================================
server = Bonito.Server(dashboard, "127.0.0.1", 8080)
println("Dashboard is running! Open http://127.0.0.1:8080 in your browser.")