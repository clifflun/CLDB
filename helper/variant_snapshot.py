import streamlit as st
import pandas as pd
import gzip
import subprocess
import streamlit.components.v1 as components

def vizCNV(chrom, start, end, refver, pr_df, pr_seg, pr_par2="NA", mo_seg="NA", fa_se="NA", \
    is_trio="FALSE", gvcf="NA", margin=25000, highlight="TRUE"):
    """
    Executes the R script 'vizCNV.R' to generate a static plot of the CNV/SV event.
    
    Args:
        chrom (str): Chromosome name (e.g., 'chr1').
        start (int): Start position.
        end (int): End position.
        refver (str): Reference genome version ('hg19' or 'hg38').
        pr_df (str): Path to proband coverage file (bed.gz).
        pr_seg (str): Path to proband segmentation file (.tsv).
        pr_par2 (str): Path to proband Parliament2 VCF (optional).
        mo_seg, fa_se (str): Maternal/Paternal segmentation files (optional).
        is_trio (str): "TRUE" or "FALSE" string flag for R.
        gvcf (str): Path to GATK VCF (optional).
        margin (int): Window size around event to display.
        highlight (str): "TRUE" to highlight the event region.
    """
    cmd = ["Rscript", 
            "vizCNV.R",
            "--chrom", f"{chrom}",
            "--from", f"{start}", 
            "--to", f"{end}", 
            "--margin", f"{margin}", 
            "--highlight", f"{highlight}",
            "--refver", f"{refver}",
            "--pr_par2", f"{pr_par2}",
            "--pr_df", f"{pr_df}", 
            "--pr_seg", f"{pr_seg}",
            "--is_trio", f"{is_trio}", 
            "--gvcf", f"{gvcf}"
            ]
    
    # Add optional arguments if they exist
    if mo_seg != "NA":
        cmd.append("--mo_seg")
        cmd.append(f"{mo_seg}")
        
    # Execute R script and raise error if it fails
    subprocess.run(cmd, check=True)


def variant_snapshot(df, dataset, pipeline, ref, db, margin):
    """
    Main function to display variant details, interactive IGV browser, and static R plots.
    
    Args:
        df (pd.DataFrame): Single-row dataframe containing variant details.
        dataset, pipeline, db: context variables for the current app state.
        ref (str): Reference genome ('hg19', 'hg38').
        margin (int): Visualization margin.
    """
    
    # --- 1. Extract Variant Details ---
    # Squeeze() converts single-element Series to scalar values
    chrom1 = df['chrom1'].squeeze()
    chrom2 = df['chrom2'].squeeze()
    pos1 = df['pos1'].squeeze()
    pos2 = df['pos2'].squeeze()
    pt_id = df['PT_ID'].squeeze()
    family = df['FAM_ID'].squeeze()
    project = df['PROJECT'].squeeze()
    SV_len = df['SV_LEN'].squeeze()
    SV_type = df['SV_TYPE'].squeeze()
    is_proband = df['IS_PROBAND'].squeeze()

    # Extract annotation columns
    L_IDR = df['L_IDR'].squeeze()
    R_IDR = df['R_IDR'].squeeze()
    L_repeatmask = df['L_repeatmask'].squeeze()
    R_repeatmask = df['R_repeatmask'].squeeze()
    
    # Normalize chromosome names for hg38
    if ref == 'hg38':
        chrom1 = 'chr' + str(chrom1)
        chrom2 = 'chr' + str(chrom2)

    # --- 2. Initialize Default Paths ---
    # Defaults prevent crashes if metadata is missing or partial
    bam_path = "NA"
    P2_path = "NA"
    MD_path = "NA"
    GATK_path = "NA"
    is_trio = "FALSE"

    # --- 3. Load Metadata based on Database Selection ---
    # Selects the correct metadata file based on the user's DB selection
    if db == 'CLDB_SR' and 'hg19' in ref: 
        meta = pd.read_csv('/CLDB/meta/meta_SR_hg19.tsv', sep='\t', index_col=False)
        meta = meta[meta['pt_id'] == pt_id]
        bam_path = meta['BAM_path'].squeeze()
        P2_path = meta['P2_path'].squeeze()
        MD_path = meta['MD_path'].squeeze()
        GATK_path = meta['GATK_path'].squeeze()
        is_trio = meta['is_trio'].squeeze()
        
    if db == 'CLDB_SR' and 'hg38' in ref: 
        meta = pd.read_csv('/CLDB/meta/meta_SR_hg38.tsv', sep='\t', index_col=False)
        meta = meta[meta['pt_id'] == pt_id]
        bam_path = meta['BAM_path'].squeeze()
        P2_path = meta['P2_path'].squeeze()
        MD_path = meta['MD_path'].squeeze()
        GATK_path = meta['GATK_path'].squeeze()
        is_trio = meta['is_trio'].squeeze()
        
    if db == 'GREGoR_SR': 
        meta = pd.read_csv('/CLDB/meta/gregor_meta.tsv', sep='\t', index_col=False)
        meta = meta[meta['pt_id'] == pt_id]
        bam_path = meta['BAM_path'].squeeze()
        P2_path = meta['P2_path'].squeeze()
        MD_path = meta['MD_path'].squeeze()
        GATK_path = meta['GATK_path'].squeeze()
        is_trio = meta['is_trio'].squeeze()
        
    if db == 'CLDB_LR':
        meta = pd.read_csv('/CLDB/meta/meta_LR.tsv', sep='\t', index_col=False)
        meta = meta[meta['pt_id'] == pt_id]
        if ref == 'hg19':
            bam_path = meta['BAM_path_hg19'].squeeze()
        elif ref == 'hg38':
            bam_path = meta['BAM_path_hg38'].squeeze()

    # --- 4. Prepare BAM/CRAM Paths for IGV ---
    # Safe handling: verify bam_path is a string (not NaN/float from empty metadata)
    if isinstance(bam_path, str):
        bam_ext = bam_path.split('.')[-1]
        if bam_ext == 'bam':
            bam_index = bam_path + '.bai'
        elif bam_ext == 'cram': 
            bam_index = bam_path + '.crai'
        else:
            bam_index = "NA"
    else:
        bam_path = "NA"
        bam_index = "NA"

    # Set repeat mask path for IGV track
    if ref == 'hg19':
        rmsk_path = '/rmsk_hg19_sorted.bed.gz'
    elif ref == 'hg38':
        rmsk_path = '/rmsk_hg38_sorted.bed.gz'

    # --- 5. Display Text Information ---
    st.write(f'BAM/CRAM source: {bam_path}')

    info1, info2 = st.columns(2)
    with info1:
        st.write(f'Patient: {pt_id}')
        st.write(f'Family: {family}')
        st.write(f'Project: {project}')
        st.write(f'Variant Type: {SV_type}')
        st.write(f'Variant length:{SV_len}')
    with info2:
        st.write(f'L_IDR: {L_IDR}')
        st.write(f'R_IDR: {R_IDR}')
        st.write(f'L_repeatmask: {L_repeatmask}')
        st.write(f'R_repeatmask: {R_repeatmask}')

    # Show all annotations in an expander
    all_annos = df.squeeze().to_dict()
    with st.expander("Full SV information"):
        col_idx = 0
        num_columns = 3
        info_cols = st.columns(num_columns)
        for i, (k,v) in enumerate(all_annos.items()):
            info_cols[col_idx].write(f'{k}: {v}')
            if (i + 1) % (len(all_annos) // num_columns + (len(all_annos) % num_columns > 0)) == 0:
                col_idx += 1
                
    st.markdown("""---""")
    st.markdown(f"<a href='https://genome.ucsc.edu/cgi-bin/hgTracks?org=human&db={ref}&position={chrom1}:{pos1}-{pos2}'>UCSC</a>", unsafe_allow_html=True)
    
    # --- 6. Render IGV Browser ---
    if bam_path != "NA":
        igv_html = f"""
            <!DOCTYPE html>
            <html lang="en">
            <head>
                <script src="https://cdn.jsdelivr.net/npm/igv@3.2.5/dist/igv.min.js"></script>
                <style>
                    #igv-container {{
                        width: 100%;
                        height: 600px;
                    }}
                </style>
            </head>
            <body>
                <div id="igv-container"></div>
                <script>
                    document.addEventListener("DOMContentLoaded", function () {{
                        var igvContainer = document.getElementById("igv-container");
                        var options = {{
                            genome: "{ref}",
                            locus: ["{chrom1}:{pos1-1000}-{pos1+1000}", "{chrom2}:{pos2-1000}-{pos2+1000}"],	
                            showCenterGuide: true,	                
                            tracks: [
                                {{
                                    name: "{pt_id}",
                                    url: "http://pnri-app17.pnri.local:3000{bam_path}",
                                    indexURL: "http://pnri-app17.pnri.local:3000{bam_index}",
                                    type: "alignment",
                                    showSoftClips: true
                                }},
                                {{
                                    name: "RepeatMask", 
                                    type: "annotation", 
                                    format: "bed", 
                                    url: "http://pnri-app17.pnri.local:3000{rmsk_path}",
                                    indexURL: "http://pnri-app17.pnri.local:3000{rmsk_path}.tbi"
                                }}
                            ]
                        }};
                        igv.createBrowser(igvContainer, options).then(function (browser) {{
                            console.log("IGV.js Ready!"); 
                        }});
                    }});
                </script>
            </body>
            </html>
            """
        components.html(igv_html, height=650)
    else:
        st.write("IGV Visualization unavailable: BAM path missing.")

    # --- 7. Generate Static CNV Plot (R Script) ---
    if 'SR' in db and SV_type != 'BND':
        if 'chr' not in str(chrom1):
            chrom = 'chr' + str(chrom1) 
        else:
            chrom = chrom1
        
        # --- ROBUST PATH PARSING ---
        # Ensure paths are valid strings before splitting to prevent "float object has no attribute split" errors
        
        # Handle P2_path
        if isinstance(P2_path, str) and ':' in P2_path:
            pr_par2 = P2_path.split(':')[1]
        else:
            pr_par2 = "NA"

        # Handle MD_path (Mosdepth)
        if isinstance(MD_path, str) and ':' in MD_path:
            pr_df = MD_path.split(':')[1]
            pr_seg = pr_df.replace('.regions.bed.gz', '_SLM_segments.tsv')
        else:
            pr_df = "NA"
            pr_seg = "NA"
        # ---------------------------
        
        if is_trio == 1:
            is_trio = "TRUE"
        else: 
            is_trio = "FALSE"

        if is_proband != 1:
            GATK_path = "NA"
        
        # Call vizCNV safely
        try:
            vizCNV(chrom, pos1, pos2, ref, pr_df, pr_seg, pr_par2=pr_par2, is_trio=is_trio, gvcf=GATK_path, margin=margin, highlight="TRUE")
            st.image(f'/tmp/r_ggplot.png', use_column_width=True)
        except subprocess.CalledProcessError as e:
            st.error(f"Visualization failed: R script returned error. {e}")