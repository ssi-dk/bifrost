window.dash_clientside = Object.assign({}, window.dash_clientside, {
    clientside: {
        large_params_function: function (largeValue1, largeValue2) {
            return someTransform(largeValue1, largeValue2);
        },
        update_selection_count: function (data) {
            return data.length
        },
        enable_selection_button: function (indices) {
            if (indices !== undefined && indices.length) {
                return false
            } else {
                return true
            }
        },
        // add_sample_to_selection: function (n_clicks, data, selected_rows, row_data) {
        //     if (selected_rows !== undefined && selected_rows.length) {
        //         selected_rows.forEach(function (e) {
        //             var row = {
        //                 "name": row_data[e].name,
        //                 "id": row_data[e].id
        //             };
        //             if (data.indexOf(row) === -1) {
        //                 data.push(row);
        //             }
        //         });
        //         return data
        //     } else {
        //         return data
        //     }
        // },
        // Passing html as output didn't work.
        // open_close_selection_modal: function (n1, n2, is_open, samples) {
        //     var tbl = document.createElement('table');
        //     for (let index = 0; index < samples.length; index++) {
        //         const element = samples[index];
        //         var tr = tbl.insertRow();
        //         var td = tr.insertCell();
        //         td.appendChild(document.createTextNode(element));
        //     }
        //     if (n1 || n2) {
        //         return [!is_open, tbl]
        //     }
        //     return [is_open, tbl]
        // }
        fill_full_virtual_name: function (name) {
            if (name === undefined) {
                return "";
            }
            name = name.toUpperCase();
            var date = new Date();
            var ye = new Intl.DateTimeFormat('en', { year: '2-digit' }).format(date)
            var mo = new Intl.DateTimeFormat('en', { month: '2-digit' }).format(date)
            var da = new Intl.DateTimeFormat('en', { day: '2-digit' }).format(date)
            return `${ye}${mo}${da}_CUSTOM_SSIVRT_${name}`
        }
    }
});