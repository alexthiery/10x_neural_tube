#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

def projectHeader() {
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_magenta = params.monochrome_logs ? '' : "\033[0;95m";

return """  -${c_dim}-----------------------------------------------------------------${c_reset}-
${c_magenta}                      __   ___ ${c_reset}
${c_magenta}                     /_ | / _ \\ ${c_reset}
${c_magenta}                      | || | | |__  __ ${c_reset}
${c_magenta}                      | || | | |\\ \\/ / ${c_reset}
${c_blue}  _______  _    ${c_reset}${c_green}  _  ${c_reset}${c_magenta} | || |_| | >  < ${c_reset}${c_blue}   ___    ___  ___    ___ ${c_reset}
${c_blue} |__   __|| |   ${c_reset}${c_green} (_)  ${c_reset}${c_magenta}|_| \\___/ /_/\\_\\ ${c_reset}${c_blue} |__ \\  / _ \\|__ \\  / _ \\ ${c_reset}
${c_blue}    | |   | |__   _   ___  _ __  _   _     ) || | | |  ) || | | | ${c_reset}
${c_blue}    | |   | '_ \\ | | / _ \\| '__|| | | |   / / | | | | / / | | | | ${c_reset}
${c_blue}    | |   | | | || ||  __/| |   | |_| |  / /_ | |_| |/ /_ | |_| | ${c_reset}
${c_blue}    |_|   |_| |_||_| \\___||_|    \\__, | |____| \\___/|____| \\___/ ${c_reset}
${c_blue}                                  __/ | ${c_reset}
${c_blue}                                 |___/ ${c_reset}
${c_magenta} alexthiery/10x_neural_tube ${c_reset}

-${c_dim}---------------------------------------------------------------${c_reset}-
    """.stripIndent()
}
