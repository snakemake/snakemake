'use strict';

class MetaData extends React.Component {
    constructor(props) {
        super(props);
    }

    render() {
        // Iterate over the metadata dictionary and add the entries to the landing page
        let metadatalist = [];
        for (const [key, value] of Object.entries(metadata)) {
            metadatalist.push(... this.innerRender(key, value));
        }
        console.log(metadatalist);
        return e(
            "ul",
            { },
            metadatalist
        )
    }

    innerRender(key, value) {
        // check if the value is an dictionary, if this is the case, recursively
        // add the elements to the landing page
        if (value.constructor == Object) {
            let result = [];
            for (const [inner_key, inner_value] of Object.entries(value)) {
                result.push(...this.innerRender(inner_key, inner_value));
            }
            return result;
        }
        return [
            e(
                ListHeading,
                { key: `${key}-heading`, text: key }
            ),
            e(
                ListItem,
                {
                    key: `${key}`,
                    className: "p-1",
                    dangerouslySetInnerHTML: { __html: value }
                }
            )
        ]
    }

}
